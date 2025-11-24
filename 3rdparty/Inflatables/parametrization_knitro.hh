#ifndef PARAMETRIZATION_KNITRO_HH
#define PARAMETRIZATION_KNITRO_HH

namespace parametrization {

template<class RParam>
void regularized_parametrization_knitro(RParam &rparam, const size_t maxIter, const std::vector<size_t> &fixedVars, double gradTol = 1e-8);

}

#endif /* end of include guard: PARAMETRIZATION_KNITRO_HH */
