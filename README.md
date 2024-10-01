# Gröbner Bases for Boolean Function Minimization

This is the corresponding tool for the journal paper "Gröbner Bases for Boolean Function Minimization" (currently under submission).

## Installation

Please install the Espresso logic minimzer. By default, on Linux, an executable `espresso` is expected in the path. On Windows, a file `Espresso.exe` is expected in the current working directory. Alternatively, you can set the `ESPRESSO_EXECUTABLE` environment variable containing a path to the Espresso executable.

Please install the latest [SageMath](https://www.sagemath.org/) version. Furthermore, please install PuLP with `python -m pip install pulp`

Lastly, install our package by typing `pip install -e .` in the folder where this README is located.

## Usage

The installed tool provides a command-line script called `groebner-min`.

It accepts [PLA files](https://user.engineering.uiowa.edu/~switchin/OldSwitching/espresso.5.html) as input. Note that, currently, only a single output is allowed. Furthermore, no don't-care inputs are supported.
A simple example of such a PLA file is the following:
```
.i 3
.o 1
000 1
010 1
011 1
110 1
111 1
.e
```

When invoked, our tool produces a small boolean formula that generates this truth table:
```sh
$ groebner-min small.pla 
formula = lambda x: ((~x[2] & ~x[0]) | x[1]) & 1
```

The final, minimized formula is printed in python syntax. The first variable in PLA corresponds to `x[0]`, the second variable is called `x[1]`, and so on.

For statistics about the run, you can use the `--stats` flag. For more detailed output, you can use the `-v` flag. 

By default, SINGULAR is used as a back end for Gröbner basis computations. You can switch to PolyBoRi by providing the `--pb` flag. Note that PolyBoRi is usally slower than SINGULAR.

For more command line flags, type `groebner-min --help`.

## Examples
Some more example PLA files can be found in the `examples/` directory.
