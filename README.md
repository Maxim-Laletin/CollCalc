# CollCalc
CollCalc is a C++ code for fast computation of collision integrals

## Make instructions

```bash
make annihilation process=Primakoff
```

"Annihilation" creates executables to compute 3-dimensional integrations for annihilation collision term (two unknown distribution functions), "co-annihilation" - 4-dimensional integrations for co-annihilation collision term (one unknown distribution functions). The name of the process file (without the extension, stored in the "process" directory) should follow `process=`. The executables are stored in "bin" directory. 

## Running 
In the existing version of the code the 'Annihilation' executable should be called with 4 arguments (in the command line), namely 'x', 'q', 'p' and the desired relative accuracy and 'Co-annihilation' - with 3 arguments ('x', 'q' and relative accuracy). If the co-annihilation code is called with 1 argument (the desired relative accuracy) the result will be the table with integrals computed for *x* and *q* values given by 'x.csv' and 'q.csv' files stored in the "parameters" directory.
For example: `./Primakoff_coann 1.0 0.5 0.1`

## Extension modules to Python
The user can incorporate collision term functions into their Python environment. To create the extension modules for a given model in the process file go to the `python` folder and run the following make instruction

```bash
make to_python process=Primakoff
```
A shared __collcalc_Primakoff.so library and the collcalc_Primakoff.py module (with the name corresponding to the specified "process") will be created in the `python/py` folder. Copy these files to the same directory as your Python project and import as follows

```python
import collcalc_Primakoff
```
The two functions available in the current version are `Annihilation(x,q,p,rel_err)` and `CoAnnihilation(x,q,rel_err)` that correspond to the executables described above.
