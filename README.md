<p align="center">
  <img src="DPC-Express_logo.png" alt="drawing" width="384" />
</p>

# DPC-Express
This is an implementation of DPC Express, which is an optimization algorithm for Density Peack Clustering (DPC) algorithm<sup>1</sup>. It accelerates DPC during tuning by removing redundant calculations while it preserves the same output as DPC. The details about DPC Express are described in the paper:

DPC Express: A Fast Drop-in Replacement of Density Peak Clustering via Redundancy Removal in Tuning

Organization of the repo:

* [DPC.cpp] - The main cpp source file for DPC Express
* [DPC.h] - The header file for DPC Express
* [main.cpp] - Sample test file for DPC Express
* [Data] - Datasets as well as codes to convert vector-based dataset to DPC-compatible input

More descriptions about the _Data_ folder please refer to the _readme_ in it.

## Compile DPC Express
```Bash
make
```

## Run DPC Express
```Bash
./DPC_Express data-file number-of-data-point number-of-clusters number-of-dc output-clustering-results: True or False output-decision-graphs: True or False

(e.g. ./DPC_Express Data/elliptical_10_2.data 500 10 16 False False)
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



git filter-branch -f --env-filter \
"GIT_AUTHOR_NAME='XA16'; GIT_AUTHOR_EMAIL='XA@XA.com'; \
GIT_COMMITTER_NAME='YMYTZSQF'; GIT_COMMITTER_EMAIL='syang16@ncsu.edu';" HEAD
