#!/usr/bin/env Rscript

# ==============================================================================
# == From "Summary" Table, launch plot generation on selected cells
# ==============================================================================

# ------------------------------------------------------------------------------
# -- Parse user input
# ------------------------------------------------------------------------------
args			= commandArgs(TRUE)
userID			= args[[1]]
analysisID		= args[[2]]
genome			= args[[3]]
bm				= args[[4]]
pseudoautosomal = args[[5]]

# --
setwd(paste('/local1/work/ginkgodev/uploads/', userID, sep=''))
maxPloidy = 6

# --
selectedCells	= read.table( paste(analysisID, '.config', sep=''), header=TRUE)
analysisType	= colnames(selectedCells)[1]

# --
raw = read.table('data', header=TRUE, sep="\t")
l = dim(raw)[1] # Number of bins
w = dim(raw)[2] # Number of samples
#
normal = sweep(raw+1, 2, colMeans(raw+1), '/')
normal2 = normal
#
GC = read.table(paste("/local1/work/ginkgodev/genomes/", genome, "/", pseudoautosomal, "/GC_", bm, sep=""), header=FALSE, sep="\t", as.is=TRUE)
#
bounds = read.table(paste("/local1/work/ginkgodev/genomes/", genome, "/", pseudoautosomal, "/bounds_", bm, sep=""), header=FALSE, sep="\t")
final  = read.table('SegCopy', header=TRUE, sep="\t")
fixed  = read.table('SegFixed', header=TRUE, sep="\t")
#
final  = final[,-c(1,2,3)]
fixed  = fixed[,-c(1,2,3)]



# --
cellIDs = c()
for(i in 1:length(selectedCells[,1]))
	cellIDs[i] = which(colnames(raw) == as.character(selectedCells[i, 1]))

if(is.null(cellIDs))
	stop("Error")


# -- Initialize color palette
cp = 3
col1 = col2 = matrix(0,3,2)
col1[1,] = c('darkmagenta', 'goldenrod')
col1[2,] = c('darkorange', 'dodgerblue')
col1[3,] = c('blue4', 'brown2')
col2[,1] = col1[,2]
col2[,2] = col1[,1]



# ------------------------------------------------------------------------------
# -- Plot Lorenz curves
# ------------------------------------------------------------------------------
if(analysisType == "lorenz")
{
	jpeg(filename=paste(analysisID, ".jpeg", sep=""), width=700, height=500)

	#
	legendNames = c("Perfect Uniformity")
	plottedFirst = 0
	for(j in 1:length(cellIDs))
	{
		k = cellIDs[j]

		nReads = sum(raw[,k])
		uniq = unique(sort(raw[,k]))
		lorenz = matrix(0, nrow=length(uniq), ncol=2)
		a = c(length(which(raw[,k]==0)), tabulate(raw[,k], nbins=max(raw[,k])))
		b = a*(0:(length(a)-1))
		for (i in 2:length(uniq))
		{
			lorenz[i,1] = sum(a[1:uniq[i]]) / l
			lorenz[i,2] = sum(b[2:uniq[i]]) / nReads
		}

		if(plottedFirst == 0)
		{
			par(mar=c(5.1, 4.1, 4.1, 18), xpd=TRUE)
			plot(lorenz, type="n", xlim=c(0,1), main="Lorenz Curve of Coverage Uniformity", xlab="Cumulative Fraction of Genome", ylab="Cumulative Fraction of Total Reads", xaxt="n", yaxt="n", cex.main=2, cex.axis=1.5, cex.lab=1.5)
		} else {
			points(lorenz, type="n")
		}


		if(plottedFirst == 0)
		{
			tu = par('usr')
			# par(xpd=FALSE)
			rect(tu[1], tu[3], tu[2], tu[4], col = "gray85")
			abline(h=seq(0,1,.1), col="white", lwd=2)
			abline(v=seq(0,1,.1), col="white", lwd=2)
			axis(side=1, at=seq(0,1,.1), tcl=.5, cex.axis=2)
			axis(side=2, at=seq(0,1,.1), tcl=.5, cex.axis=2)
			axis(side=3, at=seq(0,1,.1), tcl=.5, cex.axis=2, labels=FALSE)
			axis(side=4, at=seq(0,1,.1), tcl=.5, cex.axis=2, labels=FALSE)
			lines(c(0,1), c(0,1), lwd=2.5)
			tu = par('usr')
			# par(xpd=FALSE)
		}

		plottedFirst = 1
		try(lines(smooth.spline(lorenz), col=rainbow(length(cellIDs))[j], lwd=2.5), silent=TRUE)

		legendNames = c(legendNames, paste("Cell",selectedCells[j,1]))
	}

	legend("topright", inset=c(-0.65,0), legend=legendNames, fill=c("black", rainbow(length(cellIDs))), cex=1) #col1[cp,2]
	dev.off()
	file.create(paste(analysisID,'.done', sep=""))
}


# ------------------------------------------------------------------------------
# -- Plot GC curves
# ------------------------------------------------------------------------------
if(analysisType == "gc")
{
	jpeg(filename=paste(analysisID, ".jpeg", sep=""), width=700, height=500)

	#
	legendNames = c()
	plottedFirst = 0
	for(j in 1:length(cellIDs))
	{
		k = cellIDs[j]

		low = lowess(GC[,1], log(normal2[,k]), f=0.05)
		app = approx(low$x, low$y, GC[,1])
		cor = exp(log(normal2[,k]) - app$y)

		if(plottedFirst == 0) {
			par(mar=c(5.1, 4.1, 4.1, 18), xpd=TRUE)
			try(plot(GC[,1], log(normal2[,k]), main="GC Content vs. Bin Counts", type= "n", xlim=c(min(.3, min(GC[,1])), max(.6, max(GC[,1]))), xlab="GC content", ylab="Normalized Read Counts (log scale)", cex.main=2, cex.axis=1.5, cex.lab=1.5))
		} else {
			try(points(GC[,1], log(normal2[,k]), type="n"))
		}

		if(plottedFirst == 0)
		{
			tu = par('usr')
			rect(tu[1], tu[3], tu[2], tu[4], col = "gray85")
			abline(v=axTicks(1), col="white", lwd=2)
			abline(h=axTicks(2), col="white", lwd=2)
		}

		plottedFirst = 1
		try(points(app, col=rainbow(length(cellIDs))[j] ))

		legendNames = c(legendNames, paste("Cell",selectedCells[j,1]))
	}

	legend("topright", inset=c(-0.65,0), legend=legendNames, fill=c(rainbow(length(cellIDs))), cex=1) #col1[cp,2]
	dev.off()
	file.create(paste(analysisID,'.done', sep=""))
}


# ------------------------------------------------------------------------------
# -- Plot GC curves
# ------------------------------------------------------------------------------
if(analysisType == "cnvprofiles")
{
	library(scales)   # for alpha() opacity used in points() function

	nbCells = length(cellIDs)
	jpeg(filename=paste(analysisID, ".jpeg", sep=""), width=1000, height=200*nbCells)
	# layout(matrix(c(nbCells,1), nbCells, 1, byrow=TRUE))
	par(mfrow=c(nbCells,1)) 

	#
	rowID = 0
	for(k in cellIDs)
	{
		rowID = rowID + 1
		cat(k)
		# -- New cell
		plot(normal[,k], main=selectedCells[rowID,1], ylim=c(0, 8), type="n", xlab="Bin", ylab="Copy Number", cex.main=2, cex.axis=1.5, cex.lab=1.5)
		#
		tu = par('usr')
		par(xpd=FALSE)
		rect(tu[1], tu[3], tu[2], tu[4], col = "gray85")
		abline(h=0:19, lty=2)

		# -- Calculate CNmult (because not saved anywhere)
		CNmult = matrix(0,5,w)
		outerColsums = matrix(0, (20*(maxPloidy-1.5)+1), w)

		CNgrid = seq(1.5, maxPloidy, by=0.05)
		outerRaw = fixed[,k] %o% CNgrid
		outerRound = round(outerRaw)
		outerDiff = (outerRaw - outerRound) ^ 2
		outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
		CNmult[,k] = CNgrid[order(outerColsums[,k])[1:5]]

		# -- Plot
		flag=1
		points(normal[(0:bounds[1,2]),k]*CNmult[1,k], ylim=c(0, 6), pch=20, cex=1.5, col=alpha(col1[cp,flag], .2))
		points(final[(0:bounds[1,2]),k], ylim=c(0, 8), pch=20, cex=1.5, col=alpha(col2[cp,flag], .2))
		for (i in 1:(dim(bounds)[1]-1))
		{
			points((bounds[i,2]:bounds[(i+1),2]), normal[(bounds[i,2]:bounds[(i+1),2]),k]*CNmult[1,k], ylim=c(0, 6), pch=20, cex=1.5, col=alpha(col2[cp,flag], 0.2))
			points((bounds[i,2]:bounds[(i+1),2]), final[(bounds[i,2]:bounds[(i+1),2]),k], ylim=c(0, 8), pch=20, cex=1.5, col=alpha(col1[cp,flag], 0.2))
			if (flag == 1)
				flag = 2
			else
				flag = 1
		}
		points((bounds[(i+1),2]:l), normal[(bounds[(i+1),2]:l),k]*CNmult[1,k], ylim=c(0, 8), pch=20, cex=1.5, col=alpha(col2[cp,flag], .2))
		points((bounds[(i+1),2]:l), final[(bounds[(i+1),2]:l),k], ylim=c(0, 6), pch=20, cex=1.5, col=alpha(col1[cp,flag], .2))
	}

	dev.off()
	file.create(paste(analysisID,'.done', sep=""))
}


# ------------------------------------------------------------------------------
# -- Compare CNV profiles
# ------------------------------------------------------------------------------
if(analysisType == "cnvcompare")
{
	library(scales)   # for alpha() opacity used in points() function
        library(ggplot2)
        library(gridExtra)

        pos          = cbind(c(1,bounds[,2]), c(bounds[,2], l))

	inputNames = read.table( paste(analysisID, '.config', sep=''), header=TRUE, stringsAsFactors=FALSE)
	names = c()
	ident = c()
	ident_amp = c()
	dists_e = c()
	dists_m = c()
	spearman = c()
	pearson = c()

	#
	rowID = 1
        ref = inputNames[1,1]
        numBins = length(final[,ref])
        for(n in inputNames[2:length(inputNames[,1]),1])
	{
                rowID = rowID + 1
		startPos = unlist(gregexpr('resampled_', n, fixed=TRUE))[1]
		name = substr(n, startPos+10, nchar(n))
		endPos = unlist(gregexpr('.', name, fixed=TRUE))[1]
                #names[rowID-1] = substr(name, 1, endPos-1)
                names[rowID-1] = n

		num_match = 0
		num_match_amp = 0
		dist_euclidean = 0
		dist_manhattan = 0
		for (i in 1:numBins) {
			v1 = final[i,ref]
			v2 = final[i,n]
			if (v1 == v2)
				num_match = num_match + 1

			if ((v1 > 2 && v1 > 2) || (v1 < 2 && v2 < 2) || (v1 == 2 && v2 == 2))
				num_match_amp = num_match_amp + 1

			dist_euclidean = dist_euclidean + (v1 - v2)*(v1 - v2)
			dist_manhattan = dist_manhattan + abs(v1 - v2)
		}

		ident[rowID-1] = 100 * num_match / numBins
		ident_amp[rowID-1] = 100 * num_match_amp / numBins
		dists_e[rowID-1] = sqrt(dist_euclidean)
		dists_m[rowID-1] = dist_manhattan
		spearman[rowID-1] = cor(final[,ref], final[,n], method="spearman")
		pearson[rowID-1] = cor(final[,ref], final[,n], method="pearson")

		jpeg(filename=paste(n, "_scatter.jpeg", sep=""))
		plot(final[,ref], final[,n], type="p", xlab="CN (100%)", ylab=paste("CN (", names[rowID-1], "%)", sep=""), cex.main=1.5, cex.axis=1.5, cex.lab=1.5, col=alpha(rgb(0,0,1), 0.5), pch=19)
		dev.off()

                #Plot CN profiles
                jpeg(filename=paste(ref, "_comp_", n, ".jpeg", sep=""), width=3000, height=750)

                top=8
                rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
                rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
                cn1 = data.frame(x=which(final[,ref] != final[,n]), y=final[which(final[,ref] != final[,n]), ref])
                cn2 = data.frame(x=which(final[,ref] != final[,n]), y=final[which(final[,ref] != final[,n]), n])
                cnsame = data.frame(x=which(final[,ref] == final[,n]), y=final[which(final[,ref] == final[,n]), ref])
                #amp=data.frame(x=which(final[,k]>2), y=final[which(final[,k]>2),k])
                #del=data.frame(x=which(final[,k]<2), y=final[which(final[,k]<2),k])
                #flat=data.frame(x=which(final[,k]==2), y=final[which(final[,k]==2),k])
                anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))

                plot1 = ggplot() +
                  geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
                  geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
                  #geom_point(data=flat, aes(x=x, y=y), size=4) +
                  #geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
                  #geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
                  geom_point(data=cn1, aes(x=x, y=y), stroke=0, size=5, fill='red', color='red') +
                  geom_point(data=cn2, aes(x=x, y=y), stroke=0, size=5, fill='blue', color='blue') +
                  geom_point(data=cnsame, aes(x=x, y=y), stroke=0, size=5, fill='purple', color='purple') +
                  geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
                  scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
                  scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
                  labs(title=paste("Integer Copy Number Profiles for Sample \"", ref, "\" (red), and sample \"", n, "\" (blue)", sep=""), x="Chromosome", y="Copy Number", size=16) +
                  theme(plot.title=element_text(size=40, vjust=1.5)) +
                  theme(axis.title.x=element_text(size=40, vjust=-.05), axis.title.y=element_text(size=40, vjust=.1)) +
                  theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
                  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
                  theme(panel.background = element_rect(fill = 'gray90')) +
                  theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
                  theme(panel.grid.major.x = element_blank()) +
                  geom_vline(xintercept = c(1, l), size=.5) +
                  geom_hline(yintercept = c(-top*.1, top), size=.5)

                  grid.arrange(plot1, ncol=1)
                dev.off()

                print('')
	}

	jpeg(filename=paste("identical_bins.jpeg", sep=""))
	#plot(pcts, ident, main="Percentage of Identical Bins", type="p", ylim=c(0,100), xlab="% Reads", ylab="% Identical Bins", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
        op <- par(mar=c(10,4,4,2) + 0.1)
	barplot(ident, names.arg=names, las=2, xlab='% Reads', ylab='% Identical Bins', ylim=c(0,100))
	dev.off()

	jpeg(filename=paste("identical_amp_bins.jpeg", sep=""))
	#plot(pcts, ident_amp, main="% Bins with Matching Amp / Del", type="b", xlab="% Reads", ylab="% Identical Bins", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	barplot(ident_amp, names.arg=names, las=2, xlab='% Reads', ylab='% Same Amp/Del', ylim=c(0,100))
	dev.off()

	jpeg(filename=paste("euclidean_dist.jpeg", sep=""))
	#plot(pcts, dists_e, main="Euclidean Distance between CNVs", type="b", xlab="% Reads", ylab="Euclidean Distance", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	barplot(dists_e, names.arg=names, las=2, xlab='% Reads', ylab='Euclidean Distance')
	dev.off()

	jpeg(filename=paste("manhattan_dist.jpeg", sep=""))
	#plot(pcts, dists_m, main="Manhattan Distance between CNVs", type="b", xlab="% Reads", ylab="Manhattan Distance", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	barplot(dists_m, names.arg=names, las=2, xlab='% Reads', ylab='Manhattan Distance')
	dev.off()

	jpeg(filename=paste("spearman.jpeg", sep=""))
	#plot(pcts, spearman, main="Spearman Correlation between CNVs", type="b", ylim=c(0,1), xlab="% Reads", ylab="Spearman Correlation", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	barplot(spearman, names.arg=names, las=2, xlab='% Reads', ylab='Spearman Correlation', ylim=c(0,1))
	dev.off()

	jpeg(filename=paste("pearson.jpeg", sep=""))
	#plot(pcts, pearson, main="Pearson Correlation between CNVs", type="b", ylim=c(0,1), xlab="% Reads", ylab="Pearson Correlation", cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	barplot(pearson, names.arg=names, las=2, xlab='% Reads', ylab='Pearson Correlation', ylim=c(0,1))
	dev.off()

	#file.create(paste(analysisID,'.done', sep=""))
        cat(paste(inputNames[,1], collapse='\n'), file=paste(analysisID,'.done', sep=""))
}


# ------------------------------------------------------------------------------
# -- Plot MAD curves
# ------------------------------------------------------------------------------
if(analysisType == "mad")
{
	library(plyr)
	library(DNAcopy) #segmentation
	library(inline) #use of c++
	library(gplots) #visual plotting of tables
	library(scales)

	# Calculate MAD for selected cells
        cat("Step A\n")
	a = matrix(0, length(cellIDs), 4)
	rownames(a) <- colnames(normal[,cellIDs])
        cat("Step B\n")
	for(i in 1:length(cellIDs))
	{
		cell = cellIDs[i]
		a[i, 1] = mad(normal[-1    , cell] - normal[1:(l-1), cell])   # same as diff()
		a[i, 2] = mad(normal[-(1:2), cell] - normal[1:(l-2), cell])
		a[i, 3] = mad(normal[-(1:3), cell] - normal[1:(l-3), cell])
		a[i, 4] = mad(normal[-(1:4), cell] - normal[1:(l-4), cell])
	}
        cat("Step C\n")

	jpeg(filename=paste(analysisID, ".jpeg", sep=""), width=500, height=500)

	# Plot
	temp=cbind(a, array("", dim(a)[1]))
	mat=data.frame(off1=as.numeric(temp[,1]), off2=as.numeric(temp[,2]), off3=as.numeric(temp[,3]), off4=as.numeric(temp[,4]), ID=temp[,5])
	par(mar = c(7.0, 7.0, 5.0, 3.0))
	boxplot(mat$off1 ~ mat$ID, las=2, main="Median Absolute Deviation\nof Neighboring Bins", ylab="Median Absolute Deviation (MAD)", border=c("white"), cex.axis=1.5, cex.lab=1.5, cex.main=2, ylim=c(0, ceiling(max(a))) )

	tu <- par('usr')
	par(xpd=FALSE)
	rect(tu[1], tu[3], tu[2], tu[4], col="gray70", border="gray70", xpd=TRUE)
	rect(tu[1], tu[3], tu[2], tu[4], col="gray65", border="gray65", xpd=TRUE)
	rect(tu[1], tu[3], tu[2], tu[4], col="gray70", border="gray70", xpd=TRUE)
	rect(tu[1], tu[3], tu[2], tu[4], col=NULL)
	abline(h=seq(0,4,.5), col="white")
	par(new=TRUE)
	boxplot(mat$off1 ~ mat$ID, las=2, yaxt="n", outline=FALSE, col="#448766", names="", cex.axis=1.5, cex.lab=1.5, cex.main=2, ylim=c(0, ceiling(max(a))) )
	mtext("Median Absolute Deviation (MAD)", side=2, line=7, at=.5, cex=3)

	dev.off()
	file.create(paste(analysisID,'.done', sep=""))


}
