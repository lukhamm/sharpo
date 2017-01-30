sharpo
======

sharpo (Spherical HARmonics Projected Orbitals) is a python program for the
projection of molecular orbitals onto spherical harmonic functions at a given
center. Sharpo is written as a module for orbkit (https://orbkit.github.io/)
and thus features support for all kinds of output formats of molecular quantum
chemistry packages (GAMESS-US, GAUSSIAN, MOLPRO, Psi4, TURBOMOLE, ORCA).

Citation
--------
If you use sharpo, you use orbkit. The authors of orbkit ask you to cite their work
in the following way:

Gunter Hermann, Vincent Pohl, Jean Christophe Tremblay, Beate Paulus, Hans-Christian Hege, and Axel Schild,
"ORBKIT: A Modular Python Toolbox for Cross-Platform Postprocessing of Quantum Chemical Wavefunction Data", 
*J. Comput. Chem.* **2016**, `DOI: 10.1002/jcc.24358`__.

__ http://dx.doi.org/10.1002/jcc.24358

For more info, please visit their website.

Support
-------

If you need help or discover a bug (the code is only tested for problems and on computers in our group), please do not hesitate
to contact me.

Installation Requirements
-------------------------

To run sharpo on your system you will need to install the following Python modules:

1. Python 2.7 (should also work with later versions, http://www.python.org)
2. NumPy Library of high-level mathematical functions (http://www.numpy.org)
3. Orbkit (https://orbkit.github.io/)
   Orbkit additionally requires the following modules:
4. Cython (http://cython.org/)
5. SciPy Library of algorithms and mathematical tools (http://www.scipy.org)
6. h5py Interface to the HDF5 binary data format (http://www.h5py.org/)

Installation
------------

Once you have installed the required packeges, get a copy of sharpo either with
git or using a zip archive.

  * Using git:

    Clone the repository::

        $ git clone https://github.com/XXX

No additional installation of sharpo is necessary. Once you got your copy of the
code you can run sharpo e.g. directly from the command line in a terminal with::

$ python sharpo.py

or even on its own by typing::

$ sharpo.py

Linux 
.....

(tested on ubuntu derivatives only) If you run linux, installation is quick and simple.
You can install most of the above mentioned modules by executing the following commands::

$ apt-get install python2.7-dev python-numpy python-scipy cython python-h5py

If you want to install python3 instead of 2.7 you have to change/add the version number in the
package list.

Licence Note
------------

sharpo is free software: you can redistribute it and/or modify it under the 
terms of the GNU Lesser General Public License as published by the Free Software 
Foundation, either version 3 of the License, or any later version.

sharpo is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along 
with orbkit. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2016, Lukas Hammerschmidt.
