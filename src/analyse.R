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

biocLite("openCyto")


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
require(reshape2)
require(flowCore)
require(openCyto)
require("flowClust")
require(ggcyto)
data_dir <- "data/2018_01_08_20171212/"
data_dir <- "data/2018_04_10_20180305/"
data_dir <- "data/2018_04_10_20180320/"

fcs_files <- get_files(data_dir, ".fcs")
annotation <- read.csv(
  paste0(data_dir, "/annotation.csv"), sep = ";", header = TRUE
)
annotation <- as(annotation, "AnnotatedDataFrame")
rownames(annotation) <- fcs_files
x_raw <- read.flowSet(
  transformation = F,
  alter.names = T,
  phenoData = annotation,
  truncate_max_range = TRUE
)

# plot measurment according to time
system("mkdir -p results/gating/")
data <- apply(matrix(sampleNames(x_raw), ncol = 1), 1, FUN = function(x, fset){
    data <- data.frame(exprs(fset[[x]]))
    data <- cbind(
      x,
      1:nrow(data),
      data
    )
    return(data)
  }, fset = x_raw
)
data <- do.call(rbind, data)
names(data)[1:2] <- c("well", "step")
data.s <- data[data$step %in% sample(data$step, size = 1000),]
data.m <- melt(data.s, id.vars = c("well", "HDR.T"))
g <- ggplot(
    data.m,
    aes(x = HDR.T, y = value, group = well, color = well)
  ) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_wrap(~variable, scale = "free_y") +
  theme(legend.position="none")
print(g)
ggsave(
  filename = paste0("results/gating/time_raw.pdf"), plot = g,
  width = 29.7, height = 21, units = "cm", scale = 2
)

# gating nondebris
x_nonDebris <- x_raw
pdf("results/gating/nonDebris.pdf", width = 80.3, height = 110.7)
par(mfrow = c(16, 6))
for (i in 1:length(x)) {
  res1 <- flowClust::flowClust(x_raw[[i]], varNames=c("FSC.A", "SSC.A"), K=2, B=100)
  plot(res1, data=x_raw[[i]], level=0.95, z.cutoff=0)
  nonDebris <- which(
    getEstimates(res1)$locations[,1] == max(getEstimates(res1)$locations[,1])
  )
  Debris <- which(
    getEstimates(res1)$locations[,1] == min(getEstimates(res1)$locations[,1])
  )
  abline(v = getEstimates(res1)$locations[nonDebris,1],
         h = getEstimates(res1)$locations[nonDebris,2])
  x_nonDebris[[i]] <- split(
    x_raw[[i]],
    res1,
    population = list(
      nonDebris = nonDebris,
      Debris = Debris
    )
  )$nonDebris
}
dev.off()

# gating singlets
wf <- GatingSet(x_nonDebris)
asinhTrans <- asinhtGml2_trans()
tl <- transformerList(c("Y1.A", "B1.A"), asinhTrans)
wf <- transform(wf, tl)
add_pop(
  wf, alias = "singlets", pop = "singlets", parent = "root",
  dims = "FSC.A,FSC.H", gating_method = "singletGate",
  gating_args = "wider_gate=TRUE"
)
getPopStats(wf)

p <- ggcyto(wf, mapping = aes(x = FSC.A,y = FSC.H)) +
  geom_hex(bins = 50) +
  geom_gate("singlets") +
  ggcyto_par_set(limits = "instrument") +
  labs(title = "Singlets gate")
ggsave(
  filename = paste0("results/gating/singlets.pdf"), plot = p,
  width = 29.7, height = 21, units = "cm", scale = 2
)

# gating fluo
x_singlets <- getData(wf, "singlets")
x_fluo <- x_singlets
pdf("results/gating/fluo.pdf", width = 80.3, height = 110.7)
par(mfrow = c(16, 6))
for (i in 1:length(x)) {
  res1 <- flowClust::flowClust(x_singlets[[i]], varNames=c("Y1.A", "B1.A"), K=1, B=100)
  plot(res1, data=x_singlets[[i]], level=0.85, z.cutoff=0)
  x_fluo[[i]] <- split(
    x_singlets[[i]],
    res1,
    population = list(
      fluo = 1
    )
  )$fluo
}
dev.off()

data <- apply(matrix(sampleNames(x_fluo), ncol = 1), 1, FUN = function(x, fset){
    data <- data.frame(exprs(fset[[x]]))
    data <- cbind(
      x,
      1:nrow(data),
      data
    )
    return(data)
  }, fset = x_fluo
)
data <- do.call(rbind, data)
names(data)[1:2] <- c("well", "step")
data.s <- data[data$step %in% sample(data$step, size = 1000),]
data.m <- melt(data.s, id.vars = c("well", "HDR.T"))
g <- ggplot(
    data.m,
    aes(x = HDR.T, y = value, group = well, color = well)
  ) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_wrap(~variable, scale = "free_y") +
  theme(legend.position="none")
print(g)
ggsave(
  filename = paste0("results/gating/time_gated.pdf"), plot = g,
  width = 29.7, height = 21, units = "cm", scale = 2
)

# anova
x <- getData(wf, "nonOutliersFluo")
data <- apply(matrix(sampleNames(x), ncol = 1), 1,
  FUN = function(x, fset, infos){
    data <- data.frame(exprs(fset[[x]]))
    data <- cbind(
      x,
      1:nrow(data),
      data,
      infos[rownames(infos) %in% x, ]
    )
    return(data)
  },
  fset = x,
  infos = pData(phenoData(x))
)
data <- do.call(rbind, data)
names(data)[1:2] <- c("well", "step")
data$name <- as.factor(data$name)
data$ratio <- data$Y1.A / data$B1.A
data$drug <- relevel(data$drug, "None")
summary(data)
model <- lm(ratio ~ drug, data = data)
summary(model)
model <- lm(ratio ~ drug + Y1.A:HDR.T + B1.A:HDR.T , data = data)
summary(model)



morphGate <- norm2Filter("FSC-H", "SSC-H", filterId="MorphologyGate", scale=2)

rm_singlet <- gate_singlet(wf, ares = "FSC.A", height = "FSC.H")

smaller <- Subset(fs, morphGate)
truncateMax <- truncateTransform("truncate at 1", a=1)

fs <- getData(transform(wf, tl))
fs <- getData(wf)
p <- ggcyto(fs, aes(x = FSC.A)) 
p1 <- p + geom_histogram(bins = 60) + ggcyto_par_set(limits = "instrument")
print(p1)

p <- ggcyto(fs, aes(x = FSC.A, y = SSC.A)) +
  geom_hex(bins = 60) + 
  ggcyto_par_set(limits = "instrument")
print(p)
x11()
p <- ggcyto(fs, aes(x = Y1.A, y = B1.A)) +
  geom_hex(bins = 60) + 
  ggcyto_par_set(limits = "instrument")
print(p)
x11()
p <- ggcyto(fs, aes(x = FSC.A, y = FSC.H)) +
  geom_hex(bins = 60) + 
  ggcyto_par_set(limits = "instrument")
print(p)

require(openCyto)
gx <- GatingSet(x)
template = add_pop(
  gx, alias = "Singlets", pop = "Singlets", parent = "root",
  dims = "FSC.A,FSC.H", gating_method = "singletGate",
  gating_args = "wider_gate=TRUE,subsample_pct=0.1"
)
ggcyto(gx, mapping = aes(x = `FSC-A`,y = `FSC-H`)) +
  geom_hex(bins = 50) +
  geom_gate("Singlets") +
  xlim(c(0,2e5)) +
  labs(title = "Singlets gate")

template = rbind(
  template,
  add_pop(
    gx, alias = "nonDebris", pop = "nonDebris+",
    parent = "Singlets", dims = "FSC.A,SSC.A",
    gating_method = "boundary", collapseDataForGating = FALSE,
    gating_args = "min=40000,max=2.5e5"
  )
)
ggcyto(gx, mapping = aes(x = "FSC-A",y = "SSC-A"),
  subset = "Singlets") +
  geom_hex(bins = 100)  
  geom_gate("nonDebris")
  labs(title = "Debris gate")

require("ggplot2")
require("reshape2")

for (marker in c("Y1.A", "B1.A")) {
  for (well in levels(as.factor(data.s$well[data.s$marker %in% marker]))) {
    g <- ggplot(
        data.m[data.m$variable %in% marker & data.m$well %in% well, ],
        aes(x = step, y = value)
      ) +
      geom_smooth(se = TRUE) +
      geom_vline(xintercept = model$psi[,2] + model$psi[,3]) +
      theme_bw() +
      ggtitle(marker) +
      theme(legend.position="none")
    print(g)
  }
}

################################################################################
require("flowQ")

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
densityplot(~ FSC.A + SSC.A + Y1.A + B1.A, x)
mytransform <- transformList(c("FSC.A", "SSC.A", "Y1.A", "B1.A"), list(asinh))
densityplot(~ FSC.A + SSC.A + Y1.A + B1.A, transform(x, mytransform))
densityplot(~ FSC.A + SSC.A + Y1.A + B1.A, transform(x, mytransform), filter=lapply(c("FSC.A", "SSC.A", "Y1.A", "B1.A"), curv1Filter))

norm <- normalization(
  normFun=function(x, parameters, ...)
  warpSet(x, parameters, ...),
  parameters=c("FSC.A", "SSC.A", "Y1.A", "B1.A"),
  arguments=list(grouping="name", monwrd=TRUE),
  normalizationId="Warping"
)
normalize(transform(x, mytransform), norm)

densityplot(~ FSC.A + SSC.A + Y1.A + B1.A, transform(x, mytransform))






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
