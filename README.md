pomerol2triqs1.4
================

TRIQS wrapper around the Pomerol exact diagonalization library.  
Backport of [krivenko/pomerol2triqs](https://github.com/krivenko/pomerol2triqs) to TRIQS ver1.4.

To learn how to use this wrapper, see `example` subdir in the source directory.

Features
--------

* Diagonalization of finite fermionic models with Hamiltonians written in terms of second quantization operators.
* Calculation of single-particle Green's functions: `G(\tau)`, `G(i\omega_n)`, `G(\omega)`.
* Calculation of two-particle Green's functions: `G(\omega;\nu,\nu')` and `G(\omega;\ell,\ell')`.

Notation for the two-particle Green's functions is adopted from the
[PhD thesis of Lewin Boehnke](http://ediss.sub.uni-hamburg.de/volltexte/2015/7325/pdf/Dissertation.pdf).

Installation
------------

- Install the latest version of [Pomerol](http://aeantipov.github.io/pomerol/) exact diagonalization library.
- Install **version 1.4.2** of [TRIQS](https://triqs.ipht.cnrs.fr/1.x/install.html) library.
- `git clone https://github.com/j-otsuki/pomerol2triqs1.4.git`
- `mkdir pomerol2triqs1.4.build && cd pomerol2triqs1.4.build`
- `cmake ../pomerol2triqs1.4 -DCMAKE_BUILD_TYPE=Release -DTRIQS_PATH=<path_to_triqs_install_dir> -Dpomerol_DIR==<path_to_pomerol_install_dir>/share/pomerol`
- `make`
- `make install`

License
-------

Copyright (C) 2017 Igor Krivenko  
Copyright (C) 2018 Junya Otsuki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
