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
