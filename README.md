# CollCalc
CollCalc is a C++ code for fast computation of collision integrals

## Make instructions

`make annihilation process=Primakoff`
"Annihilation" creates executables to compute 3-dimensional integrations for annihilation collision term (two unknown distribution functions), "co-annihilation" - 4-dimensional integrations for co-annihilation collision term (one unknown distribution functions). The name of the process file (without the extension, stored in the "process" directory) should follow `process=` 

Running 'make' with 'Annihilation' or 'Co-annihilation' (see Makefile) should compile an executable file in the same directory. The specific model file (e.g. "Axion_coannihilation.cpp", should be stored in the "processes" directory) should be indicated in the corresponding line in the Makefile. 

In the existing version of the code the 'Annihilation' executable should be called with 4 arguments (in the command line), namely 'x', 'q', 'p' and the desired relative accuracy. 
For example: './Annihilation 1.0 0.1 0.2 0.1'
