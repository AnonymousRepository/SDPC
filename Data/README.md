
## Ready-for-use dataset
Due to the file size limitation in github, we only provide two small datasets (suffix is .data) here.
## Generate own datasets
We provide codes to convert a commonly used dataset format which uses a vector to represent every data point to the DPC-compatible input. Each line of the input dataset should correspond to one data point. Each dimension of the vector is separated by space. We call these input datasets raw datasets. We have provided several raw datasets (suffix is .rawdata). Users could try these raw datasets for generating DPC-compatible dataset. Also, users can use their own raw datasets.


## Compile codes
```Bash
make
```

## Run convertor
```Bash
./processBenchmark [ name of raw dataset ] [ number of dimension ]

(e.g. ./processBenchmark wine.rawData 11) 
```
default name of generated DPC-compatible dataset is [ name of raw dataset ].data
