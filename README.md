# criblejurkat R package

This repository is an R package for the automatic analaysis for fcs data for the jurkat drug crible project.

## Getting Started

These instructions will help you setup your environement to use this package
 
### Prerequisites

On ubunutu, you needd `R` and other packages :

```sh
apt-get update
apt-get install r-base libssl-dev libcurl4-openssl-dev libssh2-1-dev ca-certificates libxml2-dev r-cran-nloptr libboost-regex-dev 
```

In R, you need to install the packages `devtools`:

```R
install.packages("devtools")
```

### Installing

You can install the package directly from github with the following command:

```R
devtools::install_github("LBMC/criblejurkat")
```

### Examples

To analyse a single fcs folder (in the `set` folder):

```R
criblejurkat::set_analysis("data/set/")
```

To analyse many fcs folders (in the `sets` folder)
```R
criblejurkat::analysis("data/")
```

The analysis results will be copied in a `results/set` or `results/sets/setx` folder


## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details

