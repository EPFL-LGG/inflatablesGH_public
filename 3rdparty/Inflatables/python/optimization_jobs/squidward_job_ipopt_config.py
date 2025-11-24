from opt_config import *

name = 'isquidward'
num_threads = 1 # number of threads to use within each job (used to set {OMP,MKL}_NUM_THREADS)
#redirect_io = False

outputRootDirectory=INFLATABLES_PYROOT + '/SiggraphExamples/results_ipopt'

parameterChoices = {
    'input'                : [
                                {'name': 'squidward_fix_bdry',     'targetAttractedSheet': 'data/squidward_fix_bdry_init_highres.pkl.gz',     'fixedVars': FixedVarsBoundary    , 'uninflatedDefoInit': 'data/squidward_liftedSheetPositions_highres.txt.gz'},
                                {'name': 'squidward_fix_none',     'targetAttractedSheet': 'data/squidward_fix_none_init_highres.pkl.gz',     'fixedVars': FixedVarsNone        , 'uninflatedDefoInit': 'data/squidward_liftedSheetPositions_highres.txt.gz'},
                                {'name': 'squidward_fix_bdrywall', 'targetAttractedSheet': 'data/squidward_fix_bdrywall_init_highres.pkl.gz', 'fixedVars': FixedVarsBoundaryWall, 'uninflatedDefoInit': 'data/squidward_liftedSheetPositions_highres.txt.gz'}
                             ],
    'fittingWeight'        : [1e-3, 1e-4],
    'fusingCurveSmoothness': [FusingCurveSmoothnessParams(0, 0, 1, 1),
                              FusingCurveSmoothnessParams(0, 0, 0.5, 0.5),
                              FusingCurveSmoothnessParams(1, 0, 0.0, 0.0),
                              FusingCurveSmoothnessParams(0, 1, 0.0, 0.0)],
    'holdCPFixed'          : [True, False],
    'cbThreshold'          : [0.9,0.5]
}

cameraParamsForInput = {
    'squidward.*': {'deploy': ((-2.9991944332876064, 0.8287208791441225, 0.9943232411766028),
                               (0.3046861548028679, -0.05432817381082192, 0.9509020962232761),
                               (-0.03509330020535629, -0.09783280445341731, -0.008365205972768406)),
                     'flat': ((0.1518608344566912, -0.03246676055995215, 3.4669937818882275),
                              (0.0, 1.0, 0.0),
                              (0.1518608344566912, -0.03246676055995215, 0.0))}
}

parameterShortNames = { 'input': '', 'fittingWeight': 'fw', 'cbThreshold': 'cb', 'holdCPFixed': 'hcp', 'fusingCurveSmoothness': 'fcs'}
parameterFormatters = { 'input': lambda i: i['name'], 'fix': lambda f: f.name, 'holdCPFixed': int }
