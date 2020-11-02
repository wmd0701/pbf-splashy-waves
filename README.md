# PBS HS2020 - Project "PBF Splashy Waves"

This project by group 17 implements a position based water simulation

## Installation

### Git and CMAKE
Before we begin, you must have Git running, a distributed revision control system which you need to handin your assignments as well as keeping track of your code changes. We refer you to the online [Pro Git book](https://git-scm.com/book/en/v2) for more information. There you will also find [instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git]) on how to to install it. On Windows, we suggest using [git for windows](https://git-for-windows.github.io/).

CMake is the system this framework uses for cross-platform builds. If you are using Linux or macOS, I recommend installing it with a package manager instead of the CMake download page. E.g. on Debian/Ubuntu:
```
sudo apt-get install cmake
```
or with MacPorts on macOS:
```
sudo port install cmake.
```
On Windows, you can download it from:
[https://cmake.org/download/](https://cmake.org/download/)

### Note for linux users

Many linux distributions do not include `gcc` and the basic development tools in their default installation. On Ubuntu, you need to install the following packages:

```
sudo apt-get install build-essential libx11-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxi-dev libxmu-dev libblas-dev libxinerama-dev libxcursor-dev
```

If you are using linux with a virtual machine on Windows, it is *recommended* to use **Visual Studio** instead.

### Note for Windows users

`libigl` supports the **Microsoft Visual Studio 2015** compiler and later, in *64bit* mode. You can download *Visual Studio 2019 Community* for free from [here](https://visualstudio.microsoft.com/vs/).


### Cloning, Building and Running the Project
Before you are able to clone this project repository, you need to have an active [gitlab@ETH](https://gitlab.ethz.ch/) account.

In the next step you need to clone it to your local hard drive:
```
git clone https://gitlab.ethz.ch/pbs20-group17/pbs20.git
```
Next, cd into the newly created folder, and run the following commands inside the relevant subfolder to setup the build folder:
```
cd pbs20; mkdir build
cd build
cmake ..
```
On Windows, use the CMAKE gui with the buttons Configure and Generate.

Compile and run the executable, e.g. Ubuntu:
```
make && ./src/pbf-splashy-waves
```
Or use your favorite IDE. In case of Visual Studio, you need to open ```build/PBS.sln``` file.
Set the "pbf-splashy-waves" project as the startup project and run.