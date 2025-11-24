import sys
import gh_plotter_helpers as gh

data = [float(val) for val in sys.argv[1].split(',')]  
graph = sys.argv[2]
filename = sys.argv[3]
show_graph = True if sys.argv[4]=="True" else False
if graph == "Bending":
    gh.plot_bending_graph(data, filename, show_graph)
elif graph == "Stretching Stiffness":
    gh.plot_stretching_graph(data, filename, show_graph)