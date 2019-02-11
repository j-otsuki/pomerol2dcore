# pomerol2dcore

An interface between
[Pomerol](http://aeantipov.github.io/pomerol/) and
[DCore](https://github.com/issp-center-dev/DCore).
**pomerol2dcore** supports general local Hamiltonians with one-body and two-body terms, h0_{ij} and U_{ijkl}.
It offers an impurity solver within the Hubbard-I approximation.

## Installation

```
  $ cmake -DCMAKE_BUILD_TYPE=Release -Dpomerol_DIR=path_to_pomerol/share/pomerol/ ../pomerol2dcore
```

## How to use

Copy files in ``sample`` directory, and execute by
```
  $ pomerol2dcore params.in
```
