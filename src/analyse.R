setwd("~/projects/CribleJurkat")
devtools::load_all("src/pkg", reset = T)
data_dir <- "data/2018_01_08_20171212/"
data_dir <- "data/2018_04_10_20180305/"
data_dir <- "data/2018_04_10_20180320/"


fcs_raw <- load_data(data_dir)
fcs_nonDebris <- rm_debris(fcs_raw)
fcs_nonSinglets <- rm_nonsinglets(fcs_nonDebris)
fcs_data <- rm_nonsinglets(fcs_nonSinglets)

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
  infos = pData(phenoData(x_fluo))
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

