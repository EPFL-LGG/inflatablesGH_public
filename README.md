Inflatables Pluging for Grasshopper
===========

# Getting Started

## Building C++ Libraries
The C++ code relies on `Boost`, which must be installed separately.

The numerical solver depends on `Catamari`, which will be downloaded through the `cmake` file but require the `meson` build system. Alternatively one can choose to use the slower `CHOLMOD/UMFPACK` library, which must be installed separately.

The code also relies on several dependencies that are included as submodules:
[MeshFEM](https://github.com/MeshFEM/MeshFEM),
[libigl](https://github.com/libigl/libigl),

Finally, it includes a version of Keenan Crane's [stripe patterns code](https://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/)
modified to generate fusing curve patterns and fix a few issues with boundary handling.

### macOS
You can install all the mandatory dependencies on macOS with [MacPorts](https://www.macports.org). When installing SuiteSparse, be sure to get a version linked against `Accelerate.framework` rather than `OpenBLAS`; on MacPorts this is achieved by requesting the `accelerate` variant, which is no longer the default. Simulations will run over 2x slower under `OpenBLAS`.

```bash
# Build/version control tools, C++ code dependencies
sudo port install cmake boost ninja meson cgal4
sudo port install SuiteSparse +accelerate
# Dependencies for jupyterlab/notebooks
sudo port install python39
sudo port install npm6 nodejs16
# Dependencies for `shapely` module
sudo port install geos
```

## Obtaining and Building

Clone this repository *recursively* so that its submodules are also downloaded:

```bash
git clone --recursive git@github.com:EPFL-LGG/inflatablesGH.git
```

Build the C++ code and its Python bindings using `cmake` and your favorite
build system. For example, with [`ninja`](https://ninja-build.org):

```bash
cd Inflatables
mkdir build && cd build
cmake .. -GNinja
ninja
```
Note that Catamari's performance is affected by the build settings, so to get best performance you'll want to choose the `Release` (not `RelWithAssert`) build type and enable `MESHFEM_VECTORIZE`. 

#### Important (Temporary solution for MeshFem)
If the SuiteSparse version is newer than 5.9.0 and the following error is thrown during building: 

```bash
Error: 
const_cast from 'const value_type *' (aka 'const long *') to 'int64_t *' (aka 'long long *') is not allowed
```

Then add `#include <SuiteSparse_config.h>` in `SparseMatrices.hh`.

## Running the Jupyter Notebooks
The preferred way to interact with the inflatables code is in a Jupyter notebook,
using the Python bindings.
We recommend that you install the Python dependencies and JupyterLab itself in a
virtual environment (e.g., with [venv](https://docs.python.org/3/library/venv.html)).

```bash
pip3 install wheel # Needed if installing in a virtual environment
# Recent versions of jupyterlab and related packages cause problems:
#   JupyerLab 3.4 and later has a bug where the tab and status bar GUI
#                 remains visible after taking a viewer fullscreen
#   ipykernel > 5.5.5 clutters the notebook with stdout content
#   ipywidgets 8 and juptyerlab-widgets 3.0 break pythreejs
pip3 install jupyterlab==3.3.4 ipykernel==5.5.5 ipywidgets==7.7.2 jupyterlab-widgets==1.1.1
# If necessary, follow the instructions in the warnings to add the Python user
# bin directory (containing the 'jupyter' binary) to your PATH...

git clone https://github.com/jpanetta/pythreejs
cd pythreejs
pip3 install -e .
cd js
jupyter labextension install .

pip3 install matplotlib scipy networkx libigl setproctitle gmsh multiprocess
pip3 install shapely # dependency of the fabrication file generation
```

## Building C# Plugin
The C# plugin is compatible only with Rhino 7 and 8 for Mac with Intel processors, as well as Rhino 8 with ARM processors.

### Visual Studio for Mac (if not installed)
Download and install [Visual Studio 2022](https://visualstudio.microsoft.com/vs/mac/)

Check if `Mono` is installed with Visual Studio.
Open Visual Studio for Mac and click on 'Visual Studio' in the top menu bar. 
Select 'About Visual Studio' from the dropdown menu and, in the dialog that opens, you should see version information and installed components.
If `Mono` is not listed, download and install [Mono](https://www.mono-project.com/download/stable/)

Download the latest [RhinoVisualStudioExtensions](https://github.com/mcneel/RhinoCommonXamarinStudioAddin/releases).
Launch Visual Studio => Navigate to Visual Studio>Extensions.. => Click "Install from file" => Select the .mpack file.

Quit and Restart Visual Studio => Navigate to Extensions Studio>Add-ins..>Installed tab => Verify that RhinoCommon Plugin Support exists under the Debugging category.

## Building 
Open .sln project from `inflatablesGH/ghPlugins/InflatableSheet/` in Visual Studio and build it. This will copy all the .dll and .gha files (plugin files) in inflatablesGH/bin/isheet. The bin folder already contains the C++ library (.dylib file).

If the 'bin' folder is not referenced in Grasshopper, open Rhino, enter `GrasshopperDeveloperSettings` into the Command console, and add the path to the 'bin' folder to the Library Folders. Restart Rhino