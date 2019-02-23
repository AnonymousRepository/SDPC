
# Original DPC
This is an implementation of original DPC, which serves as the benchmark in the expereiments. This implementation is based on the Matlab codes provided by the authors of DPC<sup>1</sup>.  

Organization of the folder:

* [DPC.cpp] - The main cpp source file for DPC 
* [DPC.h] - The header file for DPC
* [main.cpp] - Sample test file for DPC


## Compile DPC 
```Bash
make
```

## Run DPC
```Bash
./DPC data-file number-of-data-point number-of-clusters number-of-dc output-clustering-results: True or False
(e.g. ./DPC ../Data/elliptical_10_2.data 500 10 8 False)
```


## Reference
[1] https://people.sissa.it/~laio/Research/Res_clustering.php
