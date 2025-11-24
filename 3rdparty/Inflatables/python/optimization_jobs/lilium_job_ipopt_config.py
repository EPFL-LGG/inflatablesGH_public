from opt_config import *

name = 'ililium'
num_threads = 1 # number of threads to use within each job (used to set {OMP,MKL}_NUM_THREADS)

outputRootDirectory=INFLATABLES_PYROOT + '/SiggraphExamples/results_ipopt'

parameterChoices = {
    'input'                : [
                                {'name': 'lilium_fix_bdry',     'targetAttractedSheet': 'data/lilium_fix_bdry_init_highres.pkl.gz',     'fixedVars': FixedVarsBoundary    , 'uninflatedDefoInit': 'data/lilium_liftedSheetPositions_highres.txt.gz'},
                                {'name': 'lilium_fix_none',     'targetAttractedSheet': 'data/lilium_fix_none_init_highres.pkl.gz',     'fixedVars': FixedVarsNone        , 'uninflatedDefoInit': 'data/lilium_liftedSheetPositions_highres.txt.gz'},
                                {'name': 'lilium_fix_bdrywall', 'targetAttractedSheet': 'data/lilium_fix_bdrywall_init_highres.pkl.gz', 'fixedVars': FixedVarsBoundaryWall, 'uninflatedDefoInit': 'data/lilium_liftedSheetPositions_highres.txt.gz'}
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
    'lilium.*': {'deploy': ((3.9524035708618674, -1.9638352398898393, 1.17750596169234),
                            (-0.3284402642919584, 0.15475658493160624, 0.9317603727419075),
                            (0.21549094439064703, 0.11721338888618847, -0.4853771087484845)),
                  'flat': ((-0.010018673275252243, -0.051311352106377314, 4.550722008110032),
                           (0.0, 1.0, 0.0),
                           (-0.010018673275252243, -0.051311352106377314, 0.0))}
}

parameterShortNames = { 'input': '', 'fittingWeight': 'fw', 'cbThreshold': 'cb', 'holdCPFixed': 'hcp', 'fusingCurveSmoothness': 'fcs'}
parameterFormatters = { 'input': lambda i: i['name'], 'fix': lambda f: f.name, 'holdCPFixed': int }
