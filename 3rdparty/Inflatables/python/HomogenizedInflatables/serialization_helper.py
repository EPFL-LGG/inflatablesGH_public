import gzip, dill 
import sys, os
# Get the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))

# Append the parent directory to the system path
sys.path.append(os.path.join(script_dir, '..'))
import inflatables_parametrization as parametrization

def save_parametrization_classes(target_surf, lscm_uv, lg, splines, default_pattern_params, num_params, rparam, path):
    dill.dump((target_surf, lscm_uv, splines, default_pattern_params, num_params, lg.getLines(), lg.alphaMin , lg.alphaMax , lg.betaMin , lg.betaMax, rparam.patternParamBounds, rparam.patternRegW, rparam.phiRegW, rparam.bendRegW, rparam.getVars()), gzip.open(path, "wb"))

def load_parametrization_classes(path):
    (target_surf, lscm_uv, splines, default_pattern_params, num_params, lines, alphaMin , alphaMax , betaMin , betaMax, patternParamBounds, patternRegW, phiRegW, bendRegW, rparam_vars) = dill.load(gzip.open(path, "rb"))

    lg = parametrization.LocalGlobalGenericParametrizer(target_surf, lscm_uv)

    lg.setLines(lines)

    lg.alphaMin = alphaMin
    lg.alphaMax = alphaMax
    lg.betaMin = betaMin
    lg.betaMax = betaMax

    rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, num_params)
    rparam.patternParamBounds = patternParamBounds
    rparam.patternRegW = patternRegW
    rparam.phiRegW = phiRegW
    rparam.bendRegW = bendRegW
    rparam.setVars(rparam_vars)
    return rparam, lg, target_surf, splines, default_pattern_params, num_params


# def save_parametrization_classes(target_surf, lscm_uv, lg, splines, default_pattern_params, num_params, rparam, pickle_path, target_surf_path):
#     if not os.path.exists(target_surf_path):
#         dill.dump(target_surf, gzip.open(target_surf_path, "wb"))
#     dill.dump((lscm_uv, splines, default_pattern_params, num_params, lg.getLines(), lg.alphaMin , lg.alphaMax , lg.betaMin , lg.betaMax, rparam.patternParamBounds, rparam.patternRegW, rparam.phiRegW, rparam.bendRegW, rparam.getVars()), gzip.open(pickle_path, "wb"))

# def load_parametrization_classes(pickle_path, target_surf_path):
#     target_surf = dill.load(gzip.open(target_surf_path, "rb"))
#     (lscm_uv, splines, default_pattern_params, num_params, lines, alphaMin , alphaMax , betaMin , betaMax, patternParamBounds, patternRegW, phiRegW, bendRegW, rparam_vars) = dill.load(gzip.open(pickle_path, "rb"))

#     lg = parametrization.LocalGlobalGenericParametrizer(target_surf, lscm_uv)

#     lg.setLines(lines)

#     lg.alphaMin = alphaMin
#     lg.alphaMax = alphaMax
#     lg.betaMin = betaMin
#     lg.betaMax = betaMax

#     rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, num_params)
#     rparam.patternParamBounds = patternParamBounds
#     rparam.patternRegW = patternRegW
#     rparam.phiRegW = phiRegW
#     rparam.bendRegW = bendRegW
#     rparam.setVars(rparam_vars)
#     return rparam, lg