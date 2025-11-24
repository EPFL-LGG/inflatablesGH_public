#include "parametrization_newton.hh"
#include <memory>

#include "regularized_parametrization.hh"
#include "generic_parametrization.hh"
#include "pattern_parametrization.hh"

namespace parametrization {

void setBoundConstraints(const RegularizedParametrizer &rparam, std::vector<NewtonProblem::BoundConstraint> &bc) {
    // Set bounds on alpha variables
    bc.reserve(rparam.numAlphaVars());
    const size_t nvar = rparam.numVars();
    const auto &vars = rparam.getVars();
    for (size_t i = rparam.alphaOffset(); i < nvar; ++i) {
        if ((vars[i] < rparam.alphaMin()) || (vars[i] > rparam.alphaMax()))
            throw std::runtime_error("Alpha bound violated " + std::to_string(i - rparam.alphaOffset()));
        bc.emplace_back(i, rparam.alphaMin(), NewtonProblem::BoundConstraint::Type::LOWER);
        bc.emplace_back(i, rparam.alphaMax(), NewtonProblem::BoundConstraint::Type::UPPER);
    }
}

void setBoundConstraints(const RegularizedGenericParametrizer &rparam, std::vector<NewtonProblem::BoundConstraint> &bc) {
    // Set bounds on alpha variables
    bc.reserve(rparam.numStretchVars());
    const size_t nvar = rparam.numVars();
    const auto &vars = rparam.getVars();
    for (size_t i = rparam.alphaOffset(); i < rparam.betaOffset(); ++i) {
        if ((vars[i] < rparam.alphaMin()) || (vars[i] > rparam.alphaMax()))
            throw std::runtime_error("Alpha bound violated " + std::to_string(i - rparam.alphaOffset()));
        bc.emplace_back(i, rparam.alphaMin(), NewtonProblem::BoundConstraint::Type::LOWER);
        bc.emplace_back(i, rparam.alphaMax(), NewtonProblem::BoundConstraint::Type::UPPER);
    }
    for (size_t i = rparam.betaOffset(); i < nvar; ++i) {
        if ((vars[i] < rparam.betaMin()) || (vars[i] > rparam.betaMax()))
            throw std::runtime_error("Beta bound violated " + std::to_string(i - rparam.betaOffset()));
        bc.emplace_back(i, rparam.betaMin(), NewtonProblem::BoundConstraint::Type::LOWER);
        bc.emplace_back(i, rparam.betaMax(), NewtonProblem::BoundConstraint::Type::UPPER);
    }
}

// For pattern parametrizer, use the patternParamBounds assigned by the user to configure the bound constraints.
void setBoundConstraints(const RegularizedPatternParametrizer &/*rparam*/, std::vector<NewtonProblem::BoundConstraint> &/*bc */) {
    // // Set bounds on pattern variables
    // const size_t nvar = rparam.numVars();
    // const size_t nt = rparam.RegularizedGenericParametrizer::mesh().numElements();
    // // Assumed numPatternVars = num_pattern_params * nt
    // if ((nvar - rparam.patternOffset())%nt != 0)
    //     throw std::runtime_error("Invalid number of pattern variables");
    // // const size_t num_pattern_vars = (nvar - rparam.patternOffset())/nt;
    // bc.reserve(nvar - rparam.patternOffset());
    // const auto &vars = rparam.getVars();
    // const size_t num_pattern_vars = (nvar - rparam.patternOffset())/nt;
    // for (size_t pi = 0; pi < num_pattern_vars; ++pi) {
    //     for (size_t ti = 0; ti < nt; ++ti) {
    //         const size_t i = rparam.patternOffset() + pi * nt + ti;
    //         if ((vars[i] < rparam.patternParamBounds()(pi, 0)) || (vars[i] > rparam.patternParamBounds()(pi, 1)))
    //             throw std::runtime_error("Pattern bound violated " + std::to_string(i - rparam.patternOffset()));
    //         bc.emplace_back(i, rparam.patternParamBounds()(pi, 0), NewtonProblem::BoundConstraint::Type::LOWER);
    //         bc.emplace_back(i, rparam.patternParamBounds()(pi, 1), NewtonProblem::BoundConstraint::Type::UPPER);
    //     }
    // }
}

void setBoundConstraints(const RegularizedParametrizerSVD &/* rparam */, std::vector<NewtonProblem::BoundConstraint> &/* bc */) {
    // No bound constraints...
}
void setBoundConstraints(const RegularizedGenericParametrizerSVD &/* rparam */, std::vector<NewtonProblem::BoundConstraint> &/* bc */) {
    // No bound constraints...
}

template<typename ParametrizationEnergy>
struct ParametrizationNewtonProblem : public NewtonProblem {
    ParametrizationNewtonProblem(ParametrizationEnergy &energy)
        : m_energy(energy), m_hessianSparsity(energy.hessianSparsityPattern()) {
            setBoundConstraints(energy, m_boundConstraints);
         }

    virtual void setVars(const Eigen::VectorXd &vars) override { m_energy.setVars(vars); }
    virtual const Eigen::VectorXd getVars() const override { return m_energy.getVars(); }
    virtual size_t numVars() const override { return m_energy.numVars(); }

    virtual Real energy() const override { return m_energy.energy(); }

    virtual Eigen::VectorXd gradient(bool /* freshIterate */ = false) const override {
        auto result = m_energy.gradient();
        return result;
    }

    virtual SuiteSparseMatrix hessianSparsityPattern() const override { /* m_hessianSparsity.fill(1.0); */ return m_hessianSparsity; }

protected:
    virtual void m_evalHessian(SuiteSparseMatrix &result, bool projectionMask) const override {
        result.setZero();
        m_energy.hessian(result, projectionMask);
    }
    virtual void m_evalMetric(SuiteSparseMatrix &result) const override {
        // TODO: mass matrix?
        result.setIdentity(true);
    }

    ParametrizationEnergy &m_energy;
    mutable SuiteSparseMatrix m_hessianSparsity;
};

template<class RParam>
ConvergenceReport regularized_parametrization_newton(RParam &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts) {
    // std::vector<NewtonProblem::BoundConstraint> bc;
    // setBoundConstraints(rparam, bc);
    auto problem = std::make_unique<ParametrizationNewtonProblem<RParam>>(rparam);
    problem->addFixedVariables(fixedVars);
    NewtonOptimizer opt(std::move(problem));
    opt.options = opts;
    return opt.optimize();
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template ConvergenceReport regularized_parametrization_newton<RegularizedGenericParametrizer   >(RegularizedGenericParametrizer    &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts);
template ConvergenceReport regularized_parametrization_newton<RegularizedParametrizer   >(RegularizedParametrizer    &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts);
template ConvergenceReport regularized_parametrization_newton<RegularizedGenericParametrizerSVD>(RegularizedGenericParametrizerSVD &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts);
template ConvergenceReport regularized_parametrization_newton<RegularizedParametrizerSVD>(RegularizedParametrizerSVD &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts);

template ConvergenceReport regularized_parametrization_newton<RegularizedPatternParametrizer   >(RegularizedPatternParametrizer  &rparam, const std::vector<size_t> &fixedVars, const NewtonOptimizerOptions &opts);

}
