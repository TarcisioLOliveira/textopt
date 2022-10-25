# textopt - optimization of cutting parameters based on fast topography simulation

Dependencies:
- [SFML](https://www.sfml-dev.org/)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

Contains code from [jdumas/mma](https://github.com/jdumas/mma).

## Installation

The code can be built using [CMake](https://cmake.org/). Cloning it into your 
computer, and provided you have the required dependencies and CMake can find them, 
it's just a matter of running:

```bash 
mkdir build 
cd build 
cmake ..  
make
make install
```

The executable is by default installed to the folder from where it was built.

Provided you have access to a terminal. Binaries are also distributed
here on GitHub via the right sidebar, under releases. The Windows binaries are
distributed with the required DLLs. For Linux, it's better to obtain them using
your distribution's package manager. This also applies for WSL.

## Usage

`textopt` can currently only be run via a command line interface, being supplied
with a configuration file, like so (for a UNIX-like shell):

```bash
./build/textopt ./examples/v0.1.0/opt.yaml
```

The YAML files contain all the necessary parameters to run the tests. There are
three kinds of tests available:

- `opt`: cutting parameter optimization
- `single`: generates a topography and related information from the cutting
parameters used as input
- `plot`: generates plot information of superficial area and roughness for a 
range of parameters. It creates a  `plot_*.txt` file inside the folder from 
where it was called, which can be plotted using the `/scripts/plot.py` script

For more information, see the example files inside the `/examples/` subfolders.
