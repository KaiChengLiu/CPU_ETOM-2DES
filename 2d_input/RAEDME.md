Files in this directory is the input files of 2DES simulation, each
file represent a conherent time and population time.

follow the instructions below:
1. modify .key-tmpl file according to your system.
2. modify gen_2d.sh, especially the "of_name"
3. run gen_2d.sh input files.


SIZE: define the system size including ground state and first excited states.

HEOM: define the HEOM configuration by (SITE_NUMBER) (TRUNCATION_LEVEL)

HAMILTONIAN: define the system Hamiltonian including ground state and first excited states.

DISORDER: define the static disorder matrix of Hamiltonian. The first two number in first row are sampling times and random seed

