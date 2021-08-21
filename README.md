# pomerol2dcore

An interface between [Pomerol](http://aeantipov.github.io/pomerol/) and [DCore](https://github.com/issp-center-dev/DCore).
**pomerol2dcore** supports general local Hamiltonians with one-body and two-body terms, h0_{ij} and U_{ijkl}.
It offers an impurity solver based on full diagonalization for arbitrary temperatures.

## Installation

```
  $ cmake -DCMAKE_BUILD_TYPE=Release -Dpomerol_DIR=/path_to_pomerol/share/pomerol/\
   -DCMAKE_INSTALL_PREFIX=/install_directory ../pomerol2dcore
  $ make
  $ make install
```
An executable file ``pomerol2dcore`` is installed in *install_directory*/bin.

## How to run

Copy files in ``sample`` directory, and execute by
```
  $ pomerol2dcore params.in
```

## Using pomerol solver in DCore

In the parameter file of DCore, specify the impurity solver as follows:
```
  [impurity_solver]
  name = pomerol
  exec_path{str} = /install_directory/bin/pomerol2dcore
  n_bath{int} = 3
  fit_gtol{float} = 1e-6
```
The case `n_bath{int} = 0` corresponds to the Hubbard-I approximation.
