require(quantmod)
require(flowStats)
require(ggcyto)
require(MASS)
require(gridExtra)
require(gplots)
source("../../src/func/functions.R")

source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("flowStats")
biocLite("flowQ")
biocLite("flowVS")
biocLite("ggcyto")


get_files <- function(path, regexp) {
  file_list <- base::list.files(
    path = path, full.names = TRUE, recursive = TRUE
  )
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(as.vector(unlist(file_list)))
}

################################################################################
setwd("~/projects/CribleJurkat")
devtools::load_all("src/pkg", reset = T)
require(flowCore)
data_dir <- "data/2018_01_08_20171212/"

fcs_files <- get_files(data_dir, ".fcs")
x <- read.flowSet(
  fcs_files,
  transformation = F,
  alter.names = T,
  phenoData = list(
    filename = "$FIL"
  )
)

require("ggplot2")
require("reshape2")
# plot measurment according to time
data <- apply(matrix(sampleNames(x), ncol = 1), 1, FUN = function(x, fset){
  data <- data.frame(exprs(fset[[x]]))
  data <- cbind(
    x,
    1:nrow(data),
    data
  )
  return(data)
}, fset = x)
data <- do.call(rbind, data)
names(data)[1:2] <- c("well", "step")
summary(data)
data.s <- data[data$step %in% sample(data$step, size = 1000),]
data.m <- melt(data.s, id.vars = c("well", "HDR.T"))

for (marker in names(data)[-c(1,2)]) {
  g <- ggplot(
      data.m[data.m$variable %in% marker, ],
      aes(x = HDR.T, y = value, group = well, color = well)
    ) +
    geom_smooth(se = FALSE) +
    theme_bw() +
    ggtitle(marker) +
    theme(legend.position="none")
    ggsave(
      filename = paste0("results/qa/time_", marker, ".pdf"), plot = g,
      width = 29.7, height = 21, units = "cm", scale = 2
    )
}
g <- ggplot(
    data.m,
    aes(x = HDR.T, y = value, group = well, color = well)
  ) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_wrap(~variable) +
  ggtitle(marker) +
  theme(legend.position="none")
ggsave(
  filename = paste0("results/qa/time_all.pdf"), plot = g,
  width = 29.7, height = 21, units = "cm", scale = 2
)


################################################################################
require("flowQ")
qa_outdir <- "results/qa/"
system(paste0("mkdir -p ", qa_outdir))

qp1 <- qaProcess.cellnumber(
  set = x,
  outdir = qa_outdir,
  two.sided = TRUE
)
qp2 <- qaProcess.marginevents(
  set = x,
  channels = c("FSC.A", "SSC.A", "Y1.A", "B1.A"),
  outdir = qa_outdir,
)
if (is.list(x)) {
  qp3 <- qaProcess.ECDFPlot(
    x,
    dyes = c("Y1.A", "B1.A"),
    outdir = qa_outdir
  )
  qp4 <- qaProcess.DensityPlot(
    x,
    dyes = c("Y1.A", "B1.A"),
    outdir = qa_outdir
  )
  qp5 <- qaProcess.2DStatsPlot(
    x,
    dyes = c("Y1.A", "B1.A"),
    outdir = qa_outdir,
    func = median
  )
  qp6 <- qaProcess.KLDistPlot(
    x,
    dyes = c("Y1.A", "B1.A"),
    outdir = qa_outdir
  )
}
x <- transform(x,
  "SSC.A"=asinh(`SSC.A`),
  "SSC.H"=asinh(`SSC.H`),
  "SSC.W"=asinh(`SSC.W`),
  "FSC.A"=asinh(`FSC.A`),
  "FSC.H"=asinh(`FSC.H`),
  "FSC.W"=asinh(`FSC.W`),
  "Y1.A"=asinh(`Y1.A`),
  "B1.A"=asinh(`B1.A`)
)
qp7 <- qaProcess.timeflow(
  set = x,
  channels = c("FSC.A", "SSC.A", "Y1.A", "B1.A"),
  outdir = qa_outdir
)
qp8 <- qaProcess.timeline(
  set = x,
  channels = c("FSC.A", "SSC.A", "Y1.A", "B1.A"),
  outdir = qa_outdir
)
url <- writeQAReport(x, processes=list(qp1, qp2, qp7, qp8), outdir=qa_outdir)
browseURL(url)


################################################################################
require(flowVS)
fcs_files <- get_files(data_dir, ".fcs")
x <- read.flowSet(
  fcs_files,
  transformation = F,
  alter.names = T
)
ggcyto(x[[1]], aes(x = FSC.H, y = FSC.W)) + geom_hex(bins = 100)
x <- transform(x,
  "SSC.A"=asinh(`SSC.A`),
  "SSC.H"=asinh(`SSC.H`),
  "SSC.W"=asinh(`SSC.W`),
  "FSC.A"=asinh(`FSC.A`),
  "FSC.H"=asinh(`FSC.H`),
  "FSC.W"=asinh(`FSC.W`),
  "Y1.A"=asinh(`Y1.A`),
  "B1.A"=asinh(`B1.A`)
)


x11()
ggcyto(x[[1]], aes(x = FSC.H, y = FSC.W)) + geom_hex(bins = 100)
x11()
fluo <- c("FSC.H", "FSC.W")
cofactors = estParamFlowVS(x, channels = fluo)
xvs = transFlowVS(x, channels = fluo, cofactors)
ggcyto(xvs[[1]], aes(x = FSC.H, y = FSC.W)) + geom_hex(bins = 100)

xvs = transFlowVS(x, channels = fluo, cofactors)
fluo <- c("Y1.A", "Y1.H", "Y1.W", "B1.A", "B1.H", "B1.W")
cofactors = estParamFlowVS(x, channels = fluo)
xvs = transFlowVS(x, channels = fluo, cofactors)

densityplot(~Y1.A+B1.A, x[[1]])
densityplot(~Y1.A+B1.A, xvs[[1]])

ggcyto(x[[1]], aes(x = FSC.H, y = FSC.W)) + geom_hex(bins = 100)

################################################################################
biocLite("flowTime")
require("flowTime")


################################################################################

library(flowViz)

library(ggcyto)
autoplot(x, "FSC.A", "SSC.A")
autoplot(x, "FSC.W", "FSC.H") + ggcyto_par_set(limits = "instrument")
autoplot(x[[1]], "SSC.W", "SSC.H") + ggcyto_par_set(limits = "instrument")
autoplot(x, "FSC.A", "B1.A")
autoplot(x, "SSC.A", "B1.A")
autoplot(x, "FSC.A", "Y1.A")
autoplot(x, "SSC.A", "Y1.A")
autoplot(x, "Y1.A", "B1.A")
autoplot(xvs, "Y1.A", "B1.A")


summary(x)
pData(phenoData(x))

library(ggcyto)
autoplot(x, "FSC-A", "SSC-A")
library(flowStats)
x <- transform(x, )
densityplot(~ FSC.A + SSC.A + V1.A + Y1.A + B1.A, x)
mytransform <- transformList(c("FSC.A", "SSC.A", "V1.A", "Y1.A", "B1.A"), list(asinh))
densityplot(~ FSC.A + SSC.A + V1.A + Y1.A + B1.A, transform(x, mytransform))
densityplot(~ FSC.A + SSC.A + V1.A + Y1.A + B1.A, transform(x, mytransform), filter=lapply(c("FSC.A", "SSC.A", "V1.A", "Y1.A", "B1.A"), curv1Filter))

norm <- normalization(
  normFun=function(x, parameters, ...)
  warpSet(x, parameters, ...),
  parameters=c("FSC.A", "SSC.A", "Y1.A", "B1.A"),
  arguments=list(grouping="name", monwrd=TRUE),
  normalizationId="Warping"
)
normalize(transform(x, mytransform), norm)

densityplot(~ FSC.A + SSC.A + V1.A + Y1.A + B1.A, transform(x, mytransform), filter=lapply(c("FSC.A", "SSC.A", "V1.A", "Y1.A", "B1.A"), curv1Filter))






load_annotation <- function(data_path) {
  annotation_path <- paste0(data_path, "annotation.csv")
  annotation <- read.table(annotation_path, h = T, sep = ";", stringsAsFactors = F)
  annotation$dapi <- as.factor(annotation$dapi)
  annotation$drug.time <- factor(paste(annotation$drug,".",  annotation$time, "UT", sep = ""))
  return(annotation)
}
annotation <- load_annotation(data_dir)

if (base::file.info(data_dir)$isdir) {
  fcs_files <- get_files(data_dir, ".fcs")
  x <- read.flowSet(fcs_files)
}
str(x)
pData(x) <- cbind(pData(x), annotation)
str(x)
