# CollCalc
CollCalc is a C++ code for fast computation of collision integrals

## Make instructions

```
make annihilation process=Primakoff
```

"Annihilation" creates executables to compute 3-dimensional integrations for annihilation collision term (two unknown distribution functions), "co-annihilation" - 4-dimensional integrations for co-annihilation collision term (one unknown distribution functions). The name of the process file (without the extension, stored in the "process" directory) should follow `process=`. The executables are stored in "bin" directory. 

## Running 
In the existing version of the code the 'Annihilation' executable should be called with 4 arguments (in the command line), namely 'x', 'q', 'p' and the desired relative accuracy and 'Co-annihilation' - with 3 arguments ('x', 'q' and relative accuracy). If the co-annihilation code is called with 1 argument (the desired relative accuracy) the result will be the table with integrals computed for *x* and *q* values given by 'x.csv' and 'q.csv' files stored in the "parameters" directory.
For example: `./Primakoff_coann 1.0 0.5 0.1`
