# Installation

On ubunutu, you needd `R` and other packages :
```sh
apt-get update
apt-get install r-base libssl-dev libcurl4-openssl-dev libssh2-1-dev ca-certificates
```

Withing a R terminal you need to install the R packages `devtools`:

```R
install.packages("devtools")
devtools::install_github("LBMC/Jalinot/criblejurkat.git")
```
