# vad
Voice activity detection based on Sohn's algorithm with the fast noise tracking by Zhang.

# Requirements
To use the library you have to install
* Boost (tested with Boost ver. 1.71)
* Armadillo (tested with Armadillo ver. 10.5.1)

# Installation
Go to root project directory and execute:

```
mkdir build
cd build
cmake ..
make
sudo make install
```

# Usage
For example usage see the main.cpp file in the example directory.
To use the library in your CMake project:

```
find_package(sohnvad REQUIRED)

add_executable( <your target> <your sources> )

target_link_libraries( <your target> sohnvad)
```
