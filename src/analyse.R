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
require("ggcyto")
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
x <- getData(x_fluo)
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

