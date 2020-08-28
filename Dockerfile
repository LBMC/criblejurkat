FROM lbmc/r-base:4.0.0

## copy files
COPY ./criblejurkat_1.2.0.tar.gz ./
RUN apk add openssl bash R-doc file \
&& Rscript -e "\
  options(repos=structure(c(CRAN=\"http://cloud.r-project.org\"))); \
  install.packages(c(\"nlme\", \"Matrix\")); \
  install.packages(\"https://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-30.tar.gz\", dep=T); \
  if (!requireNamespace(\"BiocManager\", quietly = TRUE))\
      install.packages(\"BiocManager\");\
  BiocManager::install(c('ggplot2', 'flowCore', 'flowClust', 'openCyto', 'ggcyto', 'MASS', 'flowWorkspace', 'scales', 'Biobase', 'biglm')); \
  install.packages(\"criblejurkat_1.2.0.tar.gz\", dependencies = T)"

CMD ["R"]


