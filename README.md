# EMOL

EMOL: Molecular orbital calculation program
      for educational or prototyping environment

EMOL is a small molecular orbital (MO) program written in Ruby, with
C extensions for the time-critical parts. It is meant as a prototyping
environment in which new formulations can be tried quickly.

> **Status: prototype.** EMOL is not a production code. Interfaces and
> input keys may change between versions.

## Relation to hexpr

The CSF-based Hamiltonian matrix element expressions derived from
Sasaki's spherical-tensor recoupling method were first implemented and
tested in EMOL (`src/emol_expression.rb`, `src/emol_symbolic*.rb`,
`extsrc/symbolic`). The complete formulation and a library
implementation are given in **hexpr**:

- hexpr library: <https://github.com/emol-project/hexpr>
  (DOI: [10.5281/zenodo.22974832](https://doi.org/10.5281/zenodo.22974832))
- hexpr paper (preprint): DOI [10.5281/zenodo.22977329](https://doi.org/10.5281/zenodo.22977329)

## Features

- Molecular integrals over Cartesian or spherical Gaussian basis
  functions through libcint (`extsrc/molint_libcint`)
- RHF and UHF (`src/emol_rhf.rb`, `src/emol_uhf.rb`)
- AO-to-MO integral transformation (`src/emol_trans.rb`, `extsrc/trn`)
- MP2
- Configuration interaction on a distinct row table (first- and
  second-order CI spaces; `src/emol_drt_*.rb`), with the Hamiltonian
  matrix built from numeric or symbolic matrix element expressions
  (`src/emol_hmatrix*.rb`) and an iterative eigensolver (`src/emol_liu.rb`)
- CASSCF by the super-CI method (`src/emol_superci*.rb`)
- Natural orbital and population analysis
- Basis sets (`basis/`): Sapporo DZ, DZP and TZP (2012), 6-31G(d)

## Requirements

- Ruby, with the C headers for building extensions
- C compiler
- BLAS and LAPACK (OpenBLAS preferred)
- libcint <https://github.com/sunqm/libcint> (version 3 or later)
- GNU GSL <https://www.gnu.org/software/gsl/>
- Ruby libraries: SciRuby <https://github.com/SciRuby>
  (nmatrix, narray, rb-gsl)

## Building

```sh
cd extsrc
export LIBCINT_DIR=/path/to/libcint   # if libcint is not in a standard location
./make.sh          # runs extconf.rb and make in each module
./make.sh clean
```

The integral module is built as `Molint` and copied into `src/`.
On macOS, libcint, OpenBLAS and LAPACK installed with Homebrew are found
automatically.

## Running the tutorials

```sh
cd tutorials
./calc_pop.sh      # input, integrals, RHF, population analysis
./calc_ci.sh       # ... transformation, CI, natural orbitals
./calc_superci.sh  # CASSCF by super-CI
./calc_mp2.sh
make clean_files   # remove the intermediate Data_* files
```

Each step is a separate Ruby script (`mod_*.rb`) that reads and writes
intermediate `Data_*` files in the working directory. The input is a
namelist-like format; see the `*.inp` files in `tutorials/`. In the
`$mol` section, `angular_type= 0` selects Cartesian and `angular_type= 1`
spherical Gaussian functions. Reference outputs are given for H2O and
CH4 with the Sapporo TZP basis (`*-sapTZP.out`, `*-sapTZP_sph.out`).

## Known issues

- The `extconf.rb` files of the `trn`, `wig` and `symbolic` modules
  still contain include and library paths for the original development
  machine; edit them if the build does not find GSL or the Ruby headers.
- nmatrix and rb-gsl are no longer actively maintained.

## Authors

- Hiroaki Honda
- Takeshi Noro

## License

EMOL is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version. See `COPYING.LESSER` and `COPYING`.

libcint is not included in this repository and is distributed under its
own license.

## Citation

If you use EMOL, please cite it using the metadata in `CITATION.cff`.
