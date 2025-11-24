import sys
from grasshopper_helper import generate_target_mesh
model_name = sys.argv[1]


print ('Optimizing ' + model_name)
generate_target_mesh(model_name)