<p align="center">
  <img src="SDPC_logo.png" alt="drawing" width="320" />
</p>

# DPC-Express
This is an implementation of SDPC, which is an optimization algorithm for Density Peack Clustering (DPC) algorithm<sup>1</sup>. It consistently accelerates DPC during tuning by removing redundant calculations and reduces memory consumption via stream-style data process. In the meanwhile, it preserves the same semantic of DPC and thus the same output as DPC does. The details about SDPC are described in the paper:

Streamline Density Peak Clustering for Practical Adoptions

Organization of the repo:

* [DPC.cpp] - The main cpp source file for SDPC
* [DPC.h] - The header file for SDPC
* [main.cpp] - Sample test file for SDPC
* [Data] - Datasets as well as codes to convert vector-based dataset to DPC-compatible input
* [Original] - source files for original DPC

More descriptions about the _Data_ folder and the _Original_ folder please refer to the _readme_ in them.

## Compile SDPC
```Bash
make
```

## Run SDPC
```Bash
./SDPC data-file-of-distances number-of-data-point number-of-clusters number-of-dc output-clustering-results: True or False output-decision-graphs: True or False data-file-of-data points dimension-of-data-point

(e.g. ./SDPC Data/elliptical_10_2.rawData.data 500 10 16 False False Data/elliptical_10_2.rawData 2)
```
Current version does not provide parameters to specify the range of percentage to be tested, the default setting is from 1% to 5%. Users can manually change codes in _main.cpp_ if other ranges are needed to covered. 

User can specify whether output clustering results or decision graph. In speed test, we set them as False. In practical use, they are set as True.

Clustering results of each trial will be stored in _Result_ folder and decision graphs are stored in _DecisionGraph_ folder. The format of the name of result files and decision graphs are _[trial #].res_ and _DecisionGraph.[trail #].dat_ , respectively.

Clustering results have only one column which is the cluster # of the data point, the order of data point is the same as that in the dataset.

Decision graphs are organized as two columns: 1st column is the ρ while 2nd column is the δ. It's easy to two generate graph based on this file. Users are free to use any other library or software to generate the plots. I just give an example using _gnuplot_. (Assume the decision graph file is _DecisionGraph.1.dat_)

```Bash
gnuplot -e 'plot "DecisionGraph/DecisionGraph.1.dat" with point pt 7; pause -1'
```

## Reference
[1] Rodriguez, A. and Laio, A., 2014. Clustering by fast search and find of density peaks. Science, 344(6191), pp.1492-1496.
