import sys
sys.path.append("../")
sys.path.append("../Visualization/")
sys.path.append("../../")
sys.path.append(".")
sys.path.append('../HomogenizedInflatables/');
sys.path.append('../../gmsh')
sys.path.append('parametrization_experiments')
sys.path.append('inverse_design')
import parametrization_experiments.run_parametrization_experiment as run_parametrization
import inverse_design.run_meshing_and_inflation as run_meshing_and_inflation
import inverse_design.run_inverse_design as run_inverse_design
import time 
import parallelism, multiprocessing, itertools, setproctitle

if __name__ == '__main__':
    # Change the following to a real time stamp or your own label.
    # time_stamp = time.strftime("%Y_%m_%d_%H_%M")
    # time_stamp = "2024_07_21_15_20"
    time_stamp = 'new'
    # time_stamp = 'demo_1'
    # time_stamp = 'newton'

    use_knitro = True

    args = (1, 1)
    num_thread = 8

    parallelism.set_max_num_tbb_threads(num_thread)
    parallelism.set_gradient_assembly_num_threads(num_thread)
    parallelism.set_hessian_assembly_num_threads(num_thread)

    # Coarse design
    run_parametrization.run_experiment(*args, time_stamp = time_stamp, rerun_experiment=False, use_knitro = use_knitro)
    run_meshing_and_inflation.run_experiment(*args, time_stamp = time_stamp, run_inflation = False, rerun_experiment=False, frequency = 0.2)

    # Fine-tune design
    num_iterations = 200
    run_free_boundary = False
    fix_feet = False
    run_inverse_design.run_experiment(*args, time_stamp = time_stamp, num_iterations = num_iterations, run_free_boundary = run_free_boundary, fix_feet = fix_feet)
    print("Experiment {} finished!".format(time_stamp))
