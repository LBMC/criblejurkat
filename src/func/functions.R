##############################
##### PlotChannelDensity #####
##############################

PlotChannelDensity <- function(flowset, annotation.init, files.labelling, suffix, limit = F) {
	# Plot density of a fluorescent channel for stained samples.
	#
	# Args:
	#   flowset: flowset which density is willing to be plotted.
	#   annotation.init: annotation file which contains information on samples.
	#   files.labelling: list containing files names which are stained for a given fluo.
	#   suffix: suffix within the the pdf name.
	#   limit: if TRUE, a threshold on the y axis is computed for more visibility (especially for the ratio and unstained samples).
	# Returns:
	#   pdf file with one pdf per fluorescent channel.
	tmp.files <- phenoData(flowset)$file
	annotation2 <- annotation.init[sapply(tmp.files, function(x) which(annotation.init$file == x)), ]
	for (i in seq_len(length(files.labelling))) {
		ch <- names(files.labelling)[i];
		tmp.files.labelled <- tmp.files[sapply(tmp.files, function(x) x %in% files.labelling[[ch]])]
		annotation.tmp <- annotation2[sapply(tmp.files.labelled, function(x) which(annotation2$file == x)), ]; 
		n.cl.tmp <- length(unique(annotation.tmp$cell.line)); n.fact.tmp <- levels(unique(annotation.tmp$drug.time));
		col.leg <- sapply(unique(annotation.tmp$cell.line), function(x) which(levels(factor(annotation.init$cell.line)) == x))
		col <- col.leg[annotation.tmp$cell.line]
		flowViz::flowViz.par.set("gate.density", list(lwd =  1))
		lattice::lattice.options(list("panel.error" = "warning"))
		lim <- NULL
		if (limit == T) {
			high <- 1.3*mean(fsApply(flowset[tmp.files.labelled], FUN = function(fr, chr){ dat <- exprs(fr)[, chr]; max(dat)}, chr = ch)) 
			min <- 0.8*mean(fsApply(flowset[tmp.files.labelled], FUN = function(fr, chr){ dat <- exprs(fr)[, chr]; min(dat)}, chr = ch))
			lim <- c(min, high)
		}
		d = densityplot(file~., flowset[tmp.files.labelled], xlim = lim, channels = ch, alpha = .5, fill = col, scales = list(y = list(draw = T)), col = as.numeric(annotation.tmp$drug.time), lty = as.numeric(annotation.tmp$dapi), lwd = 1, overlap = 0.2, key = list(size = 1, cex = 0.5, columns = 4, lines = list(col = c(col.leg, seq_along(n.fact.tmp), 1), type = c(rep("p", n.cl.tmp), rep("l", length(n.fact.tmp)+1)),alpha= c(rep(.5, n.cl.tmp), rep(1, length(n.fact.tmp)+1)),pch = 15,lty = c(rep(1, n.cl.tmp),rep(1,length(n.fact.tmp)),2)),text = list(c(names(col.leg), n.fact.tmp, "stained for channel"))))
		pdf(paste("./plots/", ch, suffix, ".pdf", sep = ""), height = 12, width = 5)
		print(d)
		dev.off()
	}
}


###########################
##### ScatterPlotHBin #####
###########################

ScatterPlotHBin <- function(flowset, xname, yname, ylower, yhigher, title, col.by.DAPI, dapi = dapi.channel) {
        # 2D density catter plot of yname vs xname using a binning of 64.
        #
        # Args:
        #   flowset: flowset which density is willing to be plotted.
	#   xname: name of variable to be plotted on x axis (in the phenoData of the flowset)
	#   yname: name of variable to be plotted on y axis (in the phenoData of the flowset)
	#   ylower: lower limit for y axis
	#   yhigher: upper limit for y axis
	#   title: title of the plot to save
	#   col.by.DAPI: if TRUE, color the bin according to DAPI values
	#   dapi: name of the column encoding DAPI values
        # Returns:
        #   pdf file with 2D desity scatter plot.

	if (col.by.DAPI) {
		eval(parse(text = paste("g <- ggplot(flowset, aes(x = `", xname, "`, y = `", yname, "`)) + geom_point(aes(colour = ", dapi, ")) + scale_fill_gradientn(colours = gray.colors(5))+facet_grid(cell.line~time)+geom_smooth(method = \"lm\")", sep = "")))
	} else {
		eval(parse(text = paste("g <- ggplot(flowset, aes(x = `", xname, "`, y = `", yname, "`)) + geom_hex(bins = 64) + scale_fill_gradientn(colours = gray.colors(5))+facet_grid(cell.line~time)+geom_smooth(method = \"lm\")", sep = "")))
	}
	if (is.na(ylower) == F & is.na(yhigher) == F) {
 		eval(parse(text = paste("g <- g + ylim(", ylower, ", ", yhigher, ")", sep = "")))
	}
	pdf(title, paper = "a4")
	print(g)
	dev.off()
}


#############################
##### MeanSignalsPlates #####
#############################

MeanSignalsPlates <- function(annotation.init, flowset) {
	# Compute mean signal values per channel, per well and per plate.
	#
	# Args:
	#   annotation.init: annotation file which contains information on samples.
        #   flowset: current flowset which density is willing to be plotted.
	# Returns:
	#   List of length equal to the number of plate with one value per channel, per well.
	plate.id <- unique(annotation.init$plate)
	mean.signal.plate <- list()
	for (j in plate.id) {
		ind <- which(annotation.init$plate == j); ind <- order(annotation.init$well[ind]); mean.signal <- list()
		for (i in seq_len(nb.channels)) {
			mean.signal <- c(mean.signal, list(fsApply(flowset[ind], FUN = function(fr) mean(exprs(fr)[, channels.all[i]]))))
		}
		mean.signal <- c(mean.signal, list(fsApply(flowset[ind], FUN = function(fr) {rfp <- exprs(fr)[, channels.ratio["RFP"]]; gfp <- exprs(fr)[, channels.ratio["GFP"]]; ind.rm <- c(which(rfp <= 0), which(gfp <= 0)); if(length(ind.rm)>0){res <- mean(log2(rfp[-ind.rm]/gfp[-ind.rm]))}else{res <- mean(log2(rfp/gfp))}; res})))	
		names(mean.signal) <- c(channels.all, "log2.ratio")
		mean.signal.plate <- c(mean.signal.plate, list(mean.signal))
	}
	names(mean.signal.plate) <- plate.id
	return(mean.signal.plate)
}


#################################
##### PlotMeanSignalsPlates #####
#################################

PlotMeanSignalsPlates <- function(annotation.init, name.save, mean.signal.plate) {
	# Plot mean signal values per channel, per well and per plate on pdf.
	#
	# Args:
	#   annotation.init: annotation file which contains information on samples.
	# Returns:
	#   Pdf 

	plate.id <- unique(annotation.init$plate)
	code.well <- unique(annotation.init$cpde.well)
	for (j in seq_along(plate.id)) {
		pdf(paste(name.save, j, ".pdf", sep = ""), width = 11, height= 12)
		par(mfrow = c(nb.channels+1, 2))
		tmp <- mean.signal.plate[[j]]
		for (i in seq_len(nb.channels+1)) {
			ch <- channels[i, ];temp.channels <-ch$fcs.column; if (i == nb.channels+1){temp.channels <- "log2.ratio"}
			signal <- tmp[[temp.channels]]
			info <- sapply(rownames(signal), function(x) annotation.init[annotation.init$file == x, c("well", "code.well", "cell.line", "file", "drug")])
			w <- unlist(info["well", ]); code.w <- unlist(info["code.well", ]); f <- unlist(info["file", ]); 
			t <- as.factor(unlist(info["cell.line", ])); d <- unlist(info["drug", ])
			ch.stained <- which(f %in% files.labelled[[temp.channels]])
			col <- rep("black", length(w)); col[ch.stained] <- "red"
			lwd <- rep(1, length(w)); lwd[d != "None"] <- 2.5
			ylim <- summary(tmp[[i]][, 1])[c("Min.", "Max.")]; add <- paste(" - ", temp.channels, sep = "")
			# Plot by column and by row
			for (type.plot in c("column", "row")) {
				if (type.plot == "column") {
					plot(column.template[1, code.w], tmp[[temp.channels]][, 1], col = "white", main = paste("Plate ", plate.id[j], add, " - per COLUMN", sep = ""), xaxt = "n", xlab = "Column id", ylim = ylim, xlim = c(1, 96), ylab = "mean fluo for channel")
					for (column in colnames(mat.value.template)) {
						ind <- which(unlist(sapply(code.w, function(x) length(grep(x = x, pattern = column))))!=0)
						if(length(ind)>0) {
							x.plot <- column.template[1, code.w[ind]] 
							col.tmp <- col[ind]; lwd.tmp <- lwd[ind]; pch.tmp <- as.numeric(t[ind])
							points(x.plot, tmp[[temp.channels]][, 1][ind], col = col.tmp, lwd = lwd.tmp, pch = pch.tmp, cex = 1.25)
							expr <- try(lo <- loess(tmp[[temp.channels]][, 1][ind] ~x.plot), silent = T)
							if(class(expr) != "try-error") {
							 if(length(x.plot)>2) {
							   lo <- loess(tmp[[temp.channels]][, 1][ind] ~x.plot)
							   points( x.plot, predict(lo, x.plot ),type="l", col = col.tmp[1], lty = 3)
							 }else{
							   points( x.plot, tmp[[temp.channels]][, 1][ind] ,type="l", col = col.tmp[1], lty = 3)
							 }
						  }
						}
					}
					x <- seq(8, 96, 8); 
					axis(1, at = x-4, labels = colnames(mat.value.template), lwd.ticks = 0)
					abline(v = x+0.5, lty = 2)
				}else{
					plot(row.template[1, code.w], tmp[[temp.channels]][, 1], col = "white", main = paste("Plate ", plate.id[j], add, " - per ROW", sep = ""), xaxt = "n", xlab = "Row id", ylim = ylim, xlim = c(1, 96), ylab = "mean fluo for channel")
					for (row in rownames(mat.value.template)) {
						ind <- which(unlist(sapply(unname(code.w), function(x) length(grep(x = x, pattern = row))))!=0)
						if(length(ind)>0) {
							x.plot <- row.template[1, code.w[ind]] 
							col.tmp <- col[ind]; lwd.tmp <- lwd[ind]; pch.tmp <- as.numeric(t[ind])
							points(x.plot, tmp[[temp.channels]][, 1][ind], col = col.tmp, lwd = lwd.tmp, pch = pch.tmp, cex = 1.25)
							expr <- try(lo <- loess(tmp[[temp.channels]][, 1][ind] ~x.plot), silent = T)
							if(class(try) !="try-error"){
							  if(length(x.plot)>2) {
							    lo <- loess(tmp[[temp.channels]][, 1][ind] ~x.plot)
							    points( x.plot, predict(lo, x.plot ),type="l", col = col.tmp[1], lty = 3)
							  }else{
							    points( x.plot, tmp[[temp.channels]][, 1][ind],type="l", col = col.tmp[1], lty = 3)
							  }
							}
						}
					}
					x <- seq(12, 96, 12); 
					axis(1, at = x-6, labels = rownames(mat.value.template), lwd.ticks = 0)
					abline(v = x+0.5, lty = 2)
				}
				legend("topleft", inset = c(0, -0.15), c("unlabelled", "labelled", "column (1st: A01-H01)", "loess per row"), xpd = T, horiz = T, col = c(1,2,1,2), pch = c(19, 19, -1,-1), lty = c(-1,-1,2,3), bty = "n", cex = 0.8)
				legend("topright", inset = c(0, -0.1), c("untreated", "treated"), xpd = T, horiz = T, col = 1, pch = c(5,5), lty = c(-1, -1),lwd = c(1, 2.5), bty = "n", cex = 0.8)
				legend("bottomright", inset = c(0, -0.3), xpd = T, levels(t), horiz = T, col = 1, cex = 0.75, pch = seq_along(levels(t)), bty = "n")
			}
		}
		dev.off()
	}
}
