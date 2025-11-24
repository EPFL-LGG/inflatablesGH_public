#ifndef FD_TESTS_HH
#define FD_TESTS_HH

inline void fd_tests(parametrization::RegularizedParametrizer &rparam) {
    auto origVars = rparam.getVars();
    {
        Eigen::VectorXd perturb = Eigen::VectorXd::Random(origVars.size());
        rparam.setVars(origVars + 1e-1 * perturb);
    }

    const double fd_eps = 1.e-3;
    std::cout << "RegularizedParametrizer energy:\t" << rparam.energy() << std::endl;
    std::cout << "Num flips:\t" << rparam.numFlips() << std::endl;

    auto vars = rparam.getVars();
    auto evalAt = [&](const Eigen::VectorXd &p) {
        rparam.setVars(p);
        Real e = rparam.energy();
        rparam.setVars(vars);
        return e;
    };
    Eigen::VectorXd perturb = Eigen::VectorXd::Random(vars.size());

    // perturb.head(rparam.alphaOffset()).setZero();

    Real plusEnergy  = evalAt(vars + fd_eps * perturb);
    Real minusEnergy = evalAt(vars - fd_eps * perturb);

    std::cout << "      FD delta E:\t" << (plusEnergy - minusEnergy) / (2 * fd_eps) << std::endl;
    std::cout << "Analytic delta E:\t" << rparam.gradient().dot(perturb) << std::endl;
    
    auto gradAt = [&](const Eigen::VectorXd &p) {
        rparam.setVars(p);
        auto g = rparam.gradient();
        rparam.setVars(vars);
        return g;
    };

    // auto zeroOutRotations = [&](Eigen::VectorXd &vec) { vec.tail(rparam.numVars() - rparam.numUVVars()).setZero(); };

    // perturb.tail(rparam.numVars() - rparam.numUVVars()).setZero();
    // zeroOutRotations(perturb);
    // perturb.tail(rparam.mesh().numTris()).setZero(); // clear out psi perturbation
    // perturb.segment(rparam.numVars() - 2 * rparam.mesh().numTris(), rparam.mesh().numTris()).setZero(); // clear out phi perturbation

    auto fd_delta_grad_E = ((gradAt(vars + fd_eps * perturb) - 
                             gradAt(vars - fd_eps * perturb)) / (2 * fd_eps)).eval();
    auto H = rparam.hessian();
    auto analytic_delta_grad_E = H.apply(perturb);

    // zeroOutRotations(fd_delta_grad_E);
    // zeroOutRotations(analytic_delta_grad_E);

    // fd_delta_grad_E.      tail(rparam.mesh().numTris()).setZero(); // clear out psi variables
    // analytic_delta_grad_E.tail(rparam.mesh().numTris()).setZero(); // clear out psi variables

    auto error = (analytic_delta_grad_E - fd_delta_grad_E).eval();
    Vector2D::Index loc;
    error.maxCoeff(&loc);
    std::cout << "Rel error FD delta grad E:\t" << error.norm() / fd_delta_grad_E.norm() << std::endl;
    std::cout << "Max abs error at " << loc << ": " << analytic_delta_grad_E[loc] << " vs " << fd_delta_grad_E[loc] << std::endl;

    rparam.setVars(origVars);
}

#endif /* end of include guard: FD_TESTS_HH */
