#include <MeshFEM/SparseMatrices.hh>

#include "parametrization_knitro.hh"
#include "regularized_parametrization.hh"
#include "generic_parametrization.hh"
#include "pattern_parametrization.hh"

#if HAS_KNITRO
#include <knitro.hh>
#endif // HAS_KNITRO

namespace parametrization {

#if HAS_KNITRO

void applyBoundConstraints(const RegularizedParametrizer &rparam, const Eigen::VectorXd &vars, std::vector<double> &loBounds, std::vector<double> &upBounds) {
    // Set bounds on alpha variables
    std::cout<<"RegularizedParametrizer::applyBoundConstraints"<<std::endl;
    const size_t nvar = rparam.numVars();
    for (size_t i = rparam.alphaOffset(); i < nvar; ++i) {
        if ((vars[i] < rparam.alphaMin()) || (vars[i] > rparam.alphaMax()))
            throw std::runtime_error("Alpha bound violated " + std::to_string(i - rparam.alphaOffset()));
        loBounds[i] = rparam.alphaMin();
        upBounds[i] = rparam.alphaMax();
    }
}

void applyBoundConstraints(const RegularizedGenericParametrizer &rparam, const Eigen::VectorXd &vars, std::vector<double> &loBounds, std::vector<double> &upBounds) {
    // Set bounds on alpha variables
    std::cout<<"RegularizedGenericParametrizer::applyBoundConstraints"<<std::endl;
    const size_t nvar = rparam.numVars();
    for (size_t i = rparam.alphaOffset(); i < rparam.betaOffset(); ++i) {
        if ((vars[i] < rparam.alphaMin()) || (vars[i] > rparam.alphaMax()))
            throw std::runtime_error("Alpha bound violated " + std::to_string(i - rparam.alphaOffset()));
        loBounds[i] = rparam.alphaMin();
        upBounds[i] = rparam.alphaMax();
    }
    
    // Set bounds on beta variables
    for (size_t i = rparam.betaOffset(); i < nvar; ++i) {
        if ((vars[i] < rparam.betaMin()) || (vars[i] > rparam.betaMax()))
            throw std::runtime_error("Beta bound violated " + std::to_string(i - rparam.betaOffset()));
        loBounds[i] = rparam.betaMin();
        upBounds[i] = rparam.betaMax();
    }
}

void applyBoundConstraints(const RegularizedPatternParametrizer &rparam, const Eigen::VectorXd &vars, std::vector<double> &loBounds, std::vector<double> &upBounds) {
    std::cout<<"RegularizedPatternParametrizer::applyBoundConstraints"<<std::endl;
    const size_t nvar = rparam.numVars();
    const size_t nt = rparam.RegularizedGenericParametrizer::mesh().numElements();
    // Assumed numPatternVars = num_pattern_params * nt
    if ((nvar - rparam.patternOffset())%nt != 0)
        throw std::runtime_error("Invalid number of pattern variables");
    const size_t num_pattern_vars = (nvar - rparam.patternOffset())/nt;
    for (size_t pi = 0; pi < num_pattern_vars; ++pi) {
        for (size_t ti = 0; ti < nt; ++ti) {
            const size_t i = rparam.patternOffset() + pi * nt + ti;
            if ((vars[i] < rparam.patternParamBounds()(pi, 0)) || (vars[i] > rparam.patternParamBounds()(pi, 1)))
                throw std::runtime_error("Pattern bound violated " + std::to_string(i - rparam.patternOffset()));
            loBounds[i] = rparam.patternParamBounds()(pi, 0);
            upBounds[i] = rparam.patternParamBounds()(pi, 1);
        }
    }
}

void applyBoundConstraints(const RegularizedParametrizerSVD &/* rparam */, const Eigen::VectorXd &/* vars */, std::vector<double> &/* loBounds */, std::vector<double> &/* upBounds */) {
    // No variable bounds...
}

template<typename ParametrizationEnergy>
struct KnitroParametrizationProblem : public KnitroProblem<KnitroParametrizationProblem<ParametrizationEnergy>> {
    using Base = KnitroProblem<KnitroParametrizationProblem<ParametrizationEnergy>>;

    KnitroParametrizationProblem(ParametrizationEnergy &e, const std::vector<size_t> &fixedVars)
        : Base(e.numVars(), /* num constraints */ e.numLinearInequalityConstraints()),
          energy(e),
          hessianSparsity(e.hessianSparsityPattern())
    {
        this->setObjType(KPREFIX(OBJTYPE_GENERAL));
        this->setObjGoal(KPREFIX(OBJGOAL_MINIMIZE));

        const size_t nvar = energy.numVars();

        
        // Set initial point and pin the fixed variables
        auto vars = energy.getVars();
        std::vector<Real> x_init(nvar);
        Eigen::Map<Eigen::VectorXd>(x_init.data(), x_init.size()) = vars;
        this->setXInitial(x_init);

        
        // Constraints are assumed of the form ax + by + c <= 0 where each of (a,b) will be fed through lics[ci].coeffs, (x, y)'s indices will be fed through lics[ci].vars, and c will be fed through lics[ci].constPart
        size_t nc = energy.numLinearInequalityConstraints();
        if (nc > 0) {
#if KNITRO_LEGACY
            throw std::runtime_error("Linear inequality constarints are not implemented yet for legacy Knitro");
#else
            auto lics = energy.getLinearInequalityConstraints();
            if (lics.size() != nc)
                throw std::runtime_error("Mismatch between number of linear inequality constraints and nc");
            for (size_t ci = 0; ci < nc; ++ci) {
                this->getConstraintsConstPart()  .add(ci, lics[ci].constPart);
                this->getConstraintsLinearParts().add(ci, knitro::KNLinearStructure(lics[ci].vars, lics[ci].coeffs));
            }
            this->setConUpBnds(std::vector<double>(nc, 0.0));
#endif
        }
        // ASSUMPTION: Either linear inequality constraints or bound constraints
        else {
            std::vector<double> loBounds(nvar, -KPREFIX(INFBOUND));
            std::vector<double> upBounds(nvar,  KPREFIX(INFBOUND));

            for (size_t vi : fixedVars)
                loBounds[vi] = upBounds[vi] = vars[vi];

            applyBoundConstraints(energy, vars, loBounds, upBounds);

            this->setVarLoBnds(loBounds);
            this->setVarUpBnds(upBounds);
        }

        // Inform Knitro of the Hessian sparsity pattern.
        {
            std::vector<int> hrows, hcols;
            std::vector<double> hvals;
            hessianSparsity.getIJV(hrows, hcols, hvals);
            this->setHessNnzPattern(knitro::KNSparseMatrixStructure(hrows, hcols));
        }
    }

    double evalFC(const double * x,
                        double * /* cval */,
                        double * objGrad,
                        double * /* jac */) {
        // assert(cval.size() == 0);
        // assert(jac.size() == 0);
        // assert(x.size()       == energy.numVars());
        // assert(objGrad.size() == energy.numVars());

        const size_t np = energy.numVars();
        energy.setVars(Eigen::Map<const Eigen::VectorXd>(x, np));
        double val = energy.energy();
        Eigen::Map<Eigen::VectorXd>(objGrad, np) = energy.gradient();
        return val;
    }

    // Gradient is evaluated in evaluateFC
    int evalGA(const double * /* x */, double * /* objGrad */, double * /* jac */) {
        return KPREFIX(RC_EVALFCGA);
    }

    int evalHess(const double * x, double objScalar, const double * /* lambda */,
                 double * hess) {
        // Note: Knitro gives us the Lagrange multipliers for the bound/fixed constraints in "lambda"

        assert(objScalar == 1.0);
        // assert(hess.size() == size_t(hessianSparsity.nz));
        const size_t np = energy.numVars();
        energy.setVars(Eigen::Map<const Eigen::VectorXd>(x, np));

        hessianSparsity.setZero();
        energy.hessian(hessianSparsity);
        Eigen::Map<Eigen::VectorXd>(hess, hessianSparsity.Ax.size()) = hessianSparsity.data();
        return 0;
    }

    ParametrizationEnergy &energy;
    SuiteSparseMatrix hessianSparsity;
};

template<class RParam>
void regularized_parametrization_knitro(RParam &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol) {
    KnitroParametrizationProblem<RParam> problem(rparam, fixedVars);
    
    // Create a solver - optional arguments:
    // exact first and second derivatives; no KTR_GRADOPT_* or KTR_HESSOPT_* parameter is needed.
    int hessopt = KPREFIX(HESSOPT_BFGS   ); // BFGS approximation
    hessopt = KPREFIX(HESSOPT_EXACT); // exact Hessian
    
    KnitroSolver solver(&problem, /* exact gradients */ 1, hessopt);
    // solver.setParam(KPREFIX(PARAM_HONORBNDS), KPREFIX(HONORBNDS_ALWAYS)); // always respect bounds during optimization
    solver.setParam(KPREFIX(PARAM_HONORBNDS), KPREFIX(HONORBNDS_AUTO));
    solver.setParam(KPREFIX(PARAM_BAR_MURULE), 1);
    solver.setParam(KPREFIX(PARAM_MAXIT), int(maxIter));
    solver.setParam(KPREFIX(PARAM_PRESOLVE), KPREFIX(PRESOLVE_YES));
    solver.setParam(KPREFIX(PARAM_ALGORITHM), KPREFIX(ALG_BAR_DIRECT));   // interior point with exact Hessian
    
    solver.setParam(KPREFIX(PARAM_HESSIAN_NO_F), KPREFIX(HESSIAN_NO_F_ALLOW)); // allow Knitro to call our hessvec with sigma = 0
    solver.setParam(KPREFIX(PARAM_PAR_NUMTHREADS), 12);
    
    solver.setParam(KPREFIX(PARAM_CG_MAXIT), 1000);

    solver.setParam(KPREFIX(PARAM_BAR_MAXREFACTOR), 5);
    
    // solver.setParam(KPREFIX(PARAM_ALGORITHM), KPREFIX(ALG_ACT_CG));
    // solver.setParam(KPREFIX(PARAM_ACT_QPALG), KPREFIX(ACT_QPALG_ACT_CG)); // default ended up choosing KPREFIX(ACT_QPALG_BAR_DIRECT)
    
    solver.setParam(KPREFIX(PARAM_LINSOLVER), KPREFIX(LINSOLVER_MKLPARDISO));
    solver.setParam(KPREFIX(PARAM_BAR_FEASIBLE), KPREFIX(BAR_FEASIBLE_NO));
    
    // solver.setParam(KTR_PARAM_ALGORITHM, KTR_ALG_IPDIRECT);

    solver.setParam(KPREFIX(PARAM_OPTTOL), gradTol);
    solver.setParam(KPREFIX(PARAM_OUTLEV), KPREFIX(OUTLEV_ITER_VERBOSE)); // more verbose output

    try {
        std::cout << "Running solver" << std::endl;
        int solveStatus = solver.solve();

        if (solveStatus != 0) {
            std::cout << std::endl;
            std::cout << "KNITRO failed to solve the problem, final status = ";
            std::cout << solveStatus << std::endl;
        }
    }
    catch (KnitroException &e) {
        printKnitroException(e);
        throw e;
    }
}
#else // !HAS_KNITRO

template<class RParam>
void regularized_parametrization_knitro(RParam &/* rparam */, const size_t /* maxIter */, const std::vector<size_t> &/* fixedVars */, double /* gradTol */) {
    throw std::runtime_error("Knitro wasn't found");
}

#endif // HAS_KNITRO

// Explicit function template instantiations
template void regularized_parametrization_knitro<RegularizedParametrizer   >(RegularizedParametrizer    &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol);
template void regularized_parametrization_knitro<RegularizedGenericParametrizer   >(RegularizedGenericParametrizer    &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol);
template void regularized_parametrization_knitro<RegularizedPatternParametrizer   >(RegularizedPatternParametrizer    &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol);
template void regularized_parametrization_knitro<RegularizedParametrizerSVD>(RegularizedParametrizerSVD &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol);

}
