#!/usr/bin/env Rscript

# ==============================================================================
# == Main CNV analysis script
# ==============================================================================

# ------------------------------------------------------------------------------
# -- Variables
# ------------------------------------------------------------------------------

# Config
main_dir="/local1/work/ginkgodev/scripts"

# User settings
args        = commandArgs(TRUE)
genome      = args[[1]]
user_dir    = args[[2]]
status      = args[[3]]
dat         = args[[4]]
stat        = as.numeric(args[[5]])
bm          = args[[6]]
cm          = args[[7]]
dm          = args[[8]]
cp          = as.numeric(args[[9]])
ref         = args[[10]]
f           = as.numeric(args[[11]])
facs        = args[[12]]
sex         = as.numeric(args[[13]])
bb          = as.numeric(args[[14]])
maxPloidy   = as.numeric(args[[15]])
minBinWidth = as.numeric(args[[16]])
improveBounds = as.numeric(args[[17]])

split_mode = 3
timeA <- proc.time()

outSuffix = ""
clusterCells = FALSE
if (improveBounds) {
  bm_fine     = args[[18]]
  dat_fine    = args[[19]]
  if (length(args) >= 20) {
    clusterCells = TRUE
    clustThresh = args[[20]]
    outSuffix = "_clust"
  }
} else {
  if (length(args) >= 18) {
    clusterCells = TRUE
    clustThresh = args[[18]]
    outSuffix = "_clust"
  }
}

cn_change_threshold = 0
outlier_threshold = 0.2

# Libraries
library('ctc')
library(DNAcopy) # Segmentation
library(inline)  # Use of c++
library(gplots)  # Visual plotting of tables
library(scales)
library(plyr)
library(ggplot2)
library(gridExtra)


# ------------------------------------------------------------------------------
# -- Initialize Variables & Pre-Process Data
# ------------------------------------------------------------------------------

statusFile = file( paste(user_dir, "/", status, sep="") )
writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", "<processingfile>Initializing Variables</processingfile>", "<percentdone>0</percentdone>", "<tree>clust.xml</tree>", "</status>"), statusFile)
close(statusFile)

# Load genome specific files
setwd(genome)
GC     = read.table(paste("GC_", bm, sep="")    , header=FALSE, sep="\t", as.is=TRUE)
loc    = read.table(bm                          , header=TRUE , sep="\t", as.is=TRUE)
bounds = read.table(paste("bounds_", bm, sep=""), header=FALSE, sep="\t")

if (improveBounds) {
  GC_fine     = read.table(paste("GC_", bm_fine, sep="")    , header=FALSE, sep="\t", as.is=TRUE)
  loc_fine    = read.table(bm_fine                          , header=TRUE , sep="\t", as.is=TRUE)
  bounds_fine = read.table(paste("bounds_", bm_fine, sep=""), header=FALSE, sep="\t")
}

# Load user data
setwd(user_dir)
raw    = read.table(dat, header=TRUE, sep="\t")
if (improveBounds) {
  raw_fine = read.table(dat_fine, header=TRUE, sep="\t")
  fine_mult = length(raw_fine[,1]) / length(raw[,1])
}
sink("results.txt")
cat(paste(paste(raw[,1], collapse=','), "\n", sep=','))
sink()

ploidy = rbind(c(0,0), c(0,0))
if (f == 1 | f == 2)
  ploidy = read.table(facs, header=FALSE, sep="\t", as.is=TRUE)  

if (clusterCells) {
  # Load clustering info
  clust <- readRDS(paste(user_dir, "/clust.hclust", sep=""))
  memb <- cutree(clust, h=clustThresh)
  groups <- tapply(names(memb), memb, c)
  sink(paste("groups", clustThresh, ".txt", sep=""))
  for (n in names(groups)) {
    cat(paste("group", n, '\t', paste(groups[[n]], sep=','), "\n", sep=""))
  }
  sink()

  for (n in names(groups)) {
    sink(paste("clust_group", n, ".cells", sep=""))
    cat(paste(groups[[n]], sep=','))
    sink()
  }

  # Combine raw counts in clusters
  new_raw = matrix(0,length(raw[,1]),length(groups))
  for (n in names(memb))
  {
    new_raw[,memb[n]] <- new_raw[,memb[n]] + raw[,n]
  }
  raw = new_raw
  colnames(raw) <- paste("group", names(groups), sep="")

  if (improveBounds) {
    # Combine raw counts in clusters
    new_raw_fine = matrix(0,length(raw_fine[,1]),length(groups))
    for (n in names(memb)) {
      new_raw_fine[,memb[n]] <- new_raw_fine[,memb[n]] + raw_fine[,n]
    }
    raw_fine = new_raw_fine
    colnames(raw_fine) <- paste("group", names(groups), sep="")
  }
}

# Remove bad bins
if (bb)
{
  badbins = read.table(paste(genome, "/badbins_", bm, sep=""), header=FALSE, sep="\t", as.is=TRUE)
  GC      = data.frame(GC[-badbins[,1], 1])
  loc     = loc[-badbins[,1], ]
  raw     = data.frame(raw[-badbins[,1], ])

  step  = 1
  chrom = loc[1,1]
  for (i in 1:nrow(loc))
  {
   if (loc[i,1] != chrom)
   {
     bounds[step,1] = chrom
     bounds[step,2] = i
     step           = step+1
     chrom          = loc[i,1]
    }
  }

  badbins_fine = read.table(paste(genome, "/badbins_", bm_fine, sep=""), header=FALSE, sep="\t", as.is=TRUE)
  GC_fine      = data.frame(GC_fine[-badbins_fine[,1], 1])
  loc_fine     = loc_fine[-badbins_fine[,1], ]
  raw_fine     = data.frame(raw_fine[-badbins_fine[,1], ])

  step  = 1
  chrom = loc_fine[1,1]
  for (i in 1:nrow(loc_fine))
  {
   if (loc_fine[i,1] != chrom)
   {
     bounds_fine[step,1] = chrom
     bounds_fine[step,2] = i
     step           = step+1
     chrom          = loc_fine[i,1]
    }
  }
}

# Initialize color palette
colors     = matrix(0,3,2)
colors[1,] = c('goldenrod', 'darkmagenta')
colors[2,] = c('dodgerblue', 'darkorange')
colors[3,] = c('brown2', 'blue4')

# Initialize data structures
l            = dim(raw)[1] # Number of bins
w            = dim(raw)[2] # Number of cells
n_ploidy     = length(seq(1.5, maxPloidy, by=0.05)) # Number of ploidy tests during CN inference
breaks       = matrix(0,l,w)
fixed        = matrix(0,l,w)
final        = matrix(0,l,w)
stats        = matrix(0,w,10)
if (clusterCells) {
  stats    = matrix(0,w,11)
}
CNmult       = matrix(0,n_ploidy,w)
CNerror      = matrix(0,n_ploidy,w)
outerColsums = matrix(0, (20*(maxPloidy-1.5)+1), w)
pos          = cbind(c(1,bounds[,2]), c(bounds[,2], l))

# Normalize cells
lowess.gc = function(jtkx, jtky) {
  jtklow = lowess(jtkx, log(jtky), f=0.05); 
  jtkz = approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}

normal  = sweep(raw+1, 2, colMeans(raw+1), '/')
for(k in 1:w)
{
  normal[,k] = lowess.gc( GC[,1], (raw[,k]+1)/mean(raw[,k]+1) )
}
normal2 = normal
lab     = colnames(normal)

# Prepare statistics
rownames(stats) = lab
if (clusterCells) {
  colnames(stats) = c("Cells", "Reads", "Bins", "Mean", "Var", "Disp", "Min", "25th", "Median", "75th", "Max")
} else {
  colnames(stats) = c("Reads", "Bins", "Mean", "Var", "Disp", "Min", "25th", "Median", "75th", "Max")
}

# Determine segmentation reference using dispersion (stat = 1) or reference sample (stat = 2)
if (stat == 1)
{
  F = normal[,which.min(apply(normal, 2, sd)/apply(normal,2,mean))[1]]
} else if (stat == 2) {
  R   = read.table(ref, header=TRUE, sep="\t", as.is=TRUE)
  low = lowess(GC[,1], log(R[,1]+0.001), f=0.05)
  app = approx(low$x, low$y, GC[,1])
  F   = exp(log(R[,1]) - app$y)
}

sink("timing.txt", append=TRUE)
timeB <- proc.time()
cat(paste("\n", (timeB-timeA)['elapsed'], "\n", sep=""))

# ------------------------------------------------------------------------------
# -- Process all cells
# ------------------------------------------------------------------------------

# Open output stream
sink("results.txt", append=TRUE)
cat(paste("Sample\tCopy_Number\tSoS_Predicted_Ploidy\tError_in_SoS_Approach\n", sep=""))

for(k in 1:w)
{
  # Generate basic statistics
  offset = 0
  if (clusterCells) {
    stats[k,1] = length(groups[[k]])
    offset = 1
  }
  stats[k,offset+1]  = sum(raw[,k])
  stats[k,offset+2]  = l
  stats[k,offset+3]  = round(mean(raw[,k]), digits=2)
  stats[k,offset+4]  = round(sd(raw[,k]), digits=2)
  stats[k,offset+5]  = round(stats[k,offset+4]/stats[k,offset+3], digits=2)
  stats[k,offset+6]  = min(raw[,k])
  stats[k,offset+7]  = quantile(raw[,k], c(.25))[[1]]
  stats[k,offset+8]  = median(raw[,k])
  stats[k,offset+9]  = quantile(raw[,k], c(.75))[[1]]
  stats[k,offset+10] = max(raw[,k])
}

sink("timing.txt", append=TRUE)
timeC <- proc.time()
cat(paste((timeC-timeB)['elapsed'], "\n", sep=""))
sink("results.txt", append=TRUE)

total_seg <- 0
total_plots <- 0
total_other <- 0

processSampleRefine <- function(k, loc, raw, GC, bounds, plot_sub='') {
  statusFile = file( paste(user_dir, "/", status, sep="") )
  writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>", lab[k], "</processingfile>", sep=""), paste("<percentdone>", (k*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
  close(statusFile)

  t1 <- Sys.time()

  normal = normal[,k]
  # ----------------------------------------------------------------------------
  # -- Segment data
  # ----------------------------------------------------------------------------

  # Compute log ratio between kth sample and reference
  if (stat == 0)
    lr = log2(normal)
  else
    lr = log2((normal)/(F))

  t2 <- Sys.time()
  total_other <<- total_other + (t2-t1)

  # Determine breakpoints and extract chrom/locations
  CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
  frag         = segs$output[,2:3]

  t3 <- Sys.time()
  total_seg <<- total_seg + (t3-t2)

  # Map breakpoints to kth sample
  len = dim(frag)[1]
  bps = array(0, len)
  for (j in 1:len)
    bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
  bps = sort(bps)
  bps[(len=len+1)] = l

  # Track global breakpoint locations
  seg_breaks = matrix(0,length(normal),1)
  seg_breaks[bps] = 1

  # Find outlier bins to subdivide;
  # Modify bins to contain median read count/bin within each segment
  fixed = matrix(0,length(normal),1)
  last_cn = median(normal[1:bps[2]])
  last_raw_mean = mean(raw[1:bps[2]])
  fixed[1:bps[2]] = last_cn

  outliers = c()
  breaks = c()
  j = 1
  chrom_i = 1

  zthresh = 1

  for(i in 2:(len-1)) {
    curr_cn = median(normal[bps[i]:(bps[i+1]-1)])
    curr_raw_mean = mean(raw[bps[i]:(bps[i+1]-1)])
    fixed[bps[i]:(bps[i+1]-1)] = curr_cn

    cn_change = curr_cn - last_cn

    if (loc[bps[i]-1,1] != loc[bps[i],1]) {
      next
    }

    num_b = length(bounds[,1])
    while (chrom_i <= num_b && bounds[chrom_i,2] < bps[i])
      chrom_i = chrom_i + 1

    if ((chrom_i <= num_b) && (bounds[chrom_i,2] != bps[i]) && (abs(cn_change) >= cn_change_threshold)) {
      outliers[j] = bps[i]
          bin_id = outliers[j]
          bin_start = loc[bin_id-2,2]
          bin_end = loc[bin_id,2]

          fine_start = 1
          found_chrom = FALSE
          for (chrom_id in 1:length(bounds_fine[,1])) {
            if (bounds_fine[chrom_id,1] == loc[bin_id,1]) {
              found_chrom = TRUE
              if (chrom_id == 1)
                fine_start = 1
              else
                fine_start = bounds_fine[chrom_id-1,2]
              break
            }
          }
          if (!found_chrom) {
            cat(paste("Couldn't find chromosome ", loc[bin_id,1], "\n", sep=""))
            quit()
          }
          while (loc_fine[fine_start,1] != loc[bin_id,1] || loc_fine[fine_start,2] <= bin_start) # Find the right bin
            fine_start = fine_start+1

          fine_end = fine_start
          while (loc_fine[fine_end+1,2] < bin_end)
            fine_end = fine_end+1

          left_bin_overlap = (loc_fine[fine_start,2] - bin_start) / (loc_fine[fine_start,2] - loc_fine[fine_start-1,2])
          right_bin_overlap = (bin_end - loc_fine[fine_end,2]) / (loc_fine[fine_end+1,2] - loc_fine[fine_end,2])

          total_bins = fine_end - fine_start
          total_bins_olap = total_bins + left_bin_overlap + right_bin_overlap
          left_bins = 0
          right_bins = total_bins
          left_count = 0
          right_count = sum(raw_fine[(fine_start+1):fine_end,k])

          min_score = (right_count*fine_mult/right_bins - curr_raw_mean) ** 2
          peak_left_bins = left_bins

          margin = 0
          for (left_bins in 1:(total_bins-margin)) {
            right_bins = total_bins - left_bins
            left_count = left_count + raw_fine[fine_start+left_bins,k]
            right_count = right_count - raw_fine[fine_start+left_bins,k]
            zscore1 = (left_count*fine_mult/left_bins - last_raw_mean) ** 2
            curr_score = zscore1
            if (right_bins > 0) {
                zscore2 = (right_count*fine_mult/right_bins - curr_raw_mean) ** 2
                curr_score = (zscore1 + zscore2)/2
            }

            if (margin == left_bins || curr_score < min_score) {
                min_score = curr_score
                peak_left_bins = left_bins
            }
          }
          new_bound = loc_fine[fine_start+peak_left_bins,2]
          loc[bin_id-1,2] = new_bound
      j = j+1
    }

    last_cn = curr_cn
    last_raw_mean = curr_raw_mean
  }
  normal = normal/mean(fixed)
  fixed = fixed/mean(fixed)

  # ----------------------------------------------------------------------------
  # -- Determine Copy Number (SoS Method)
  # ----------------------------------------------------------------------------

  # Determine Copy Number     
  CNgrid           = seq(1.5, maxPloidy, by=0.05)
  outerRaw         = fixed %o% CNgrid
  outerRound       = round(outerRaw)
  outerDiff        = (outerRaw - outerRound) ^ 2
  outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult[,k]       = CNgrid[order(outerColsums[,k])]
  CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

  if (f == 0 | length(which(lab[k]==ploidy[,1]))==0 ) {
    CN = CNmult[1,k]
  } else if (f == 1) {
    CN = ploidy[which(lab[k]==ploidy[,1]),2]
  } else {
    estimate = ploidy[which(lab[k]==ploidy[,1]),2]
    CN = CNmult[which(abs(CNmult[,k] - estimate)<.4),k][1]
  }
  final = round(fixed*CN)

  # Output results of CN calculations to file
  out=paste(lab[k], CN, paste(CNmult[,k], collapse= ","), paste(CNerror[,k], collapse= ","), sep="\t")
  cat(out, "\n")

  write.table(cbind(loc,final,raw), file=paste(lab[k], "_CN", plot_sub, ".tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc),lab[k], "Raw"), sep="\t", quote=FALSE)

  t4 <- Sys.time()
  total_other <<- total_other + (t4-t3)

  # ----------------------------------------------------------------------------
  # -- Generate Plots & Figures
  # ----------------------------------------------------------------------------

  # Plot Distribution of Read Coverage
  jpeg(filename=paste(lab[k], "_dist", plot_sub, ".jpeg", sep=""), width=3000, height=750)
  
  top=round(quantile(raw, c(.995))[[1]])
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  main=data.frame(x=which(raw<top), y=raw[which(raw<top)])
  dist_outliers=data.frame(x=which(raw>top), y=array(top*.99, length(which(raw>top))))
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))

  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) + 
    geom_point(data=main, aes(x=x, y=y), size=3) +
    geom_point(data=dist_outliers, aes(x=x, y=y), shape=5, size=4) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    labs(title=paste("Genome Wide Read Distribution for Sample \"", lab[k], "\"", sep=""), x="Chromosome", y="Read Count", size=16) +
    theme(plot.title=element_text(size=36, vjust=1.5)) +
    theme(axis.title.x=element_text(size=40, vjust=-.1), axis.title.y=element_text(size=40, vjust=-.06)) +
    theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'gray90')) +
    theme(plot.margin=unit(c(.5,1,.5,1.5),"cm")) +
    theme(panel.grid.major.x = element_blank()) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    geom_vline(xintercept = c(1, l), size=.25) +
    geom_hline(yintercept = c(-top*.1, top), size=.25)

    grid.arrange(plot1, ncol=1)
  dev.off()

  #Plot histogram of bin counts
  jpeg(filename=paste(lab[k], "_counts", plot_sub, ".jpeg", sep=""), width=2500, height=1500)
    par(mar = c(7.0, 7.0, 7.0, 3.0))

    temp=sort(raw)[round(l*.01) : (l-round(l*.01))] 
    reads = hist(temp, breaks=100, plot=FALSE)
    # plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", xaxt="n", cex.main=3, cex.axis=2, cex.lab=2)
    plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", cex.main=3, cex.axis=2, cex.lab=2)
    # axis(side=1, at=seq(min(temp), round(diff(range(temp))/20)*22, round(diff(range(temp))/20)), cex.axis=2)
    tu = par('usr')
    par(xpd=FALSE)
    clip(tu[1], as.integer(mean(temp)-(diff(reads$mids)/2)), tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(mean(temp)+(diff(reads$mids)/2), tu[2], tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(tu[1], mean(temp) - sd(temp), tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(mean(temp) + sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(tu[1], mean(temp) - 2*sd(temp), tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE)
    clip(mean(temp) + 2*sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE) 
    legend("topright", inset=.05, legend=c("mean", "< 1σ", "> 1σ", "> 2σ"), fill=c("black", "gray50", "gray75", "gray90"), cex=2.5)
  dev.off()

  #Plot lorenz curves
  jpeg(filename=paste(lab[k], "_lorenz", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  nReads=sum(raw)
  uniq=unique(sort(raw))
  
  lorenz=matrix(0, nrow=length(uniq), ncol=2)
  a=c(length(which(raw==0)), tabulate(raw, nbins=max(raw)))
  b=a*(0:(length(a)-1))
  for (i in 2:length(uniq)) {
    lorenz[i,1]=sum(a[1:uniq[i]])/l
    lorenz[i,2]=sum(b[2:uniq[i]])/nReads
  }

  # smooth.spline needs >= 4 points...
  fit = data.frame(x=lorenz[,1], y=lorenz[,2])
  if(nrow(lorenz) >= 4)
  {
    spline = try(smooth.spline(lorenz))
    if(class(spline) != "try-error")
      fit = data.frame(x=spline$x, y=spline$y)
  }

  perf=data.frame(x=c(0,1), y=c(0,1))

  plot1 = try(ggplot() +
    geom_line(data=perf, aes(x=x, y=y, color="Perfect Uniformity"), size=3) +
    geom_line(data=fit, aes(x=x, y=y, color="Sample Uniformity"), size=3) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    labs(title=paste("Lorenz Curve of Coverage Uniformity for Sample ", lab[k], sep=""), x="Cumulative Fraction of Genome", y="Cumulative Fraction of Total Reads") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) + theme(legend.title=element_blank(), legend.text=element_text(size=40)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.15, .85)) +
    scale_color_manual(name='', values=c('Perfect Uniformity'="black", 'Sample Uniformity'=colors[cp,1])))

    grid.arrange(plot1, ncol=1)
  dev.off()

  ##Plot GC correction
  #jpeg(filename=paste(lab[k], "_GC", outSuffix, ".jpeg", sep=""), width=2500, height=1250)

  #low = lowess(GC[,1], log(normal2[,k]), f=0.05)
  #app = approx(low$x, low$y, GC[,1])
  #cor = exp(log(normal2[,k]) - app$y)
  
  #uncorrected = data.frame(x=GC[,1], y=log(normal2[,k]))
  #corrected = data.frame(x=GC[,1], y=log(cor))
  #fit = data.frame(x=app$x, y=app$y)

  #try(plot1 <- ggplot() +
  #  geom_point(data=uncorrected, aes(x=x, y=y), size=3) +
  #  geom_line(data=fit, aes(x=x, y=y, color="Lowess Fit"), size=3) +
  #  scale_x_continuous(limits=c(min(.3, min(GC[,1])), max(.6, max(GC[,1]))), breaks=seq(.3,.6,.05)) +
  #  labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Uncorrected)", sep=""), x="GC content", y="Normalized Read Counts (log scale)") +
  #  theme(plot.title=element_text(size=45, vjust=1.5)) +
  #  theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
  #  theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
  #  theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
  #  theme(panel.background = element_rect(color = 'black')) +
  #  theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
  #  theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
  #  theme(legend.position=c(.85, .9)) +
  #  scale_color_manual(name='', values=colors[cp,1]))

  #try(plot2 <- ggplot() +
  #  geom_point(data=corrected, aes(x=x, y=y), size=3) +
  #  scale_x_continuous(limits=c(min(.3, min(GC[,1])), max(.6, max(GC[,1]))), breaks=seq(.3,.6,.05)) +
  #  labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Corrected)", sep=""), x="GC content", y="") +
  #  theme(plot.title=element_text(size=45, vjust=1.5)) +
  #  theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
  #  theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
  #  theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
  #  theme(panel.background = element_rect(color = 'black')))

  #  try(grid.arrange(plot1, plot2, ncol=2))
  #dev.off()

  #Plot Scaled/Normalized Bin Count Histogram
  jpeg(filename=paste(lab[k], "_hist", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  clouds=data.frame(x=normal*CN)

  plot1 = ggplot() +
    geom_histogram(data=clouds, aes(x=x), binwidth=.05, color="black", fill="gray60") +
    geom_vline(xintercept=seq(0,10,1), size=1, linetype="dashed", color=colors[cp,1]) +
    scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1)) +
    labs(title=paste("Frequency of Bin Counts for Sample \"", lab[k], "\"\nNormalized and Scaled by Predicted CN (", CNmult[1,k], ")", sep=""), x="Copy Number", y="Frequency") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black'))

    grid.arrange(plot1, ncol=1)
  dev.off()

  #Plot sum of squares error for each potential copy number
  jpeg(filename=paste(lab[k], "_SoS", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  top = max(outerColsums[,k])
  dat = data.frame(x=CNgrid, y=outerColsums[,k])
  lim = cbind(c(seq(0,5000,500), 1000000), c(50, 100, 100, 200, 250, 400, 500, 500, 600, 600, 750, 1000))
  step = lim[which(top<lim[,1])[1],]
  minSoS = data.frame(x=CNmult[1,k], y=CNerror[1,k])
  bestSoS = data.frame(x=CN, y=outerColsums[which(CNgrid==CN),k])

  plot1 = ggplot() +
    geom_line(data=dat, aes(x=x, y=y), size=3) +
    geom_point(data=dat, aes(x=x, y=y), shape=21, fill="black", size=5) +
    geom_point(data=minSoS, aes(x=x, y=y*1.02, color="Minimum SoS Error"), shape=18, size=15) +
    geom_point(data=bestSoS, aes(x=x, y=y*.98, color="Chosen Ploidy"), shape=18, size=15) +
    scale_x_continuous(limits=c(1.5, 6), breaks=seq(1.5, 6, .5)) +
    scale_y_continuous(limits=c(.5*min(outerColsums[,k]), top), breaks=seq(0, step[1], step[2])) +
    labs(title="Sum of Squares Error Across Potential Copy Number States", x="Copy Number Multiplier", y="Sum of Squares Error") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.85, .9)) +
    scale_color_manual(name='', values=c('Minimum SoS Error'=colors[cp,1], 'Chosen Ploidy'=colors[cp,2]))

    grid.arrange(plot1, ncol=1)
  dev.off()

  #Plot colored CN profile
  jpeg(filename=paste(lab[k], "_CN", plot_sub, ".jpeg", sep=""), width=3000, height=750)

  top=8
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  clouds=data.frame(x=1:length(normal), y=normal*CN)
  amp=data.frame(x=which(final>2), y=final[which(final>2)])
  del=data.frame(x=which(final<2), y=final[which(final<2)])
  flat=data.frame(x=which(final==2), y=final[which(final==2)])
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))
  outlier_pts=data.frame(x=outliers, y=normal[outliers]*CN)

  dat = data.frame(x=bps, y=final[bps])

  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
    geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
    geom_point(data=outlier_pts, aes(x=x, y=y), size=4, color='blue4') +
    #geom_point(data=flat, aes(x=x, y=y), size=4) +
    #geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
    #geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
    geom_line(data=dat, aes(x=x, y=y), size=3) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n Predicted Ploidy = ", CN, sep=""), x="Chromosome", y="Copy Number", size=16) +
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

  t5 <- Sys.time()
  total_plots <<- total_plots + (t5-t4)

  v <- list("out" = outliers, "breaks" = breaks, "loc" = loc, "gc" = GC, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
  return(v)
}
#
processSample <- function(k, loc, raw, GC, bounds, plot_sub=NULL, outliers=NULL, breaks=NULL, bin_wgts=NULL) {
  statusFile = file( paste(user_dir, "/", status, sep="") )
  writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>", lab[k], "</processingfile>", sep=""), paste("<percentdone>", (k*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
  close(statusFile)

  t1 <- Sys.time()

  cat(paste('Processing', '\n', sep=''))

  # ----------------------------------------------------------------------------
  # -- Split outlier bins
  # ----------------------------------------------------------------------------
  normal = normal[,k]

  if (!missing(outliers)) {
    num_o = length(outliers)
    if (num_o > 0) {

      if (missing(bin_wgts) || is.null(bin_wgts)) {
        bin_wgts = c()
        for (i in 1:length(raw))
          bin_wgts[i] = 1
      }

      for (i in 1:num_o) {
        bin_id = outliers[i]+i-1
        if (loc[bin_id-1,1] != loc[bin_id,1])
          next 

        bin_start = loc[bin_id-1,2]
        bin_end = loc[bin_id,2]
        bin_mid = (bin_end - bin_start) * breaks[i] + bin_start

        fine_start = 1
        found_chrom = FALSE
        for (chrom_id in 1:length(bounds_fine[,1])) {
          if (bounds_fine[chrom_id,1] == loc[bin_id,1]) {
            found_chrom = TRUE
            if (chrom_id == 1)
              fine_start = 1
            else
              fine_start = bounds_fine[chrom_id-1,2]
            break
          }
        }
        if (!found_chrom) {
          cat(paste("Couldn't find chromosome ", loc[bin_id,1], "\n", sep=""))
          quit()
        }
        #while (loc_fine[fine_start,1] != loc[bin_id+i-1,1]) # Find the right chromosome
        #  fine_start = fine_start+1
        while (loc_fine[fine_start,1] != loc[bin_id,1] || loc_fine[fine_start+1,2] <= bin_start) # Find the right bin
          fine_start = fine_start+1

        fine_mid = fine_start
        while (loc_fine[fine_mid,2] < bin_mid) {
          fine_mid = fine_mid+1
        }
        if ((bin_mid - loc_fine[fine_mid-1,2]) < (loc_fine[fine_mid,2] - bin_mid))
          fine_mid = fine_mid - 1
        #fine_mid_first = fine_mid

        fine_end = fine_start
        while (loc_fine[fine_end,2] < bin_end)
          fine_end = fine_end+1

        if ((fine_end-fine_start+1) %% 2 == 1)
          fine_mid = (fine_end + fine_start) / 2
        else {
          fine_midA = (fine_end + fine_start - 1) / 2
          fine_midB = fine_midA + 1
          if ((bin_mid - fine_midA) < (fine_midB - bin_mid))
            fine_mid = fine_midA
          else
            fine_mid = fine_midB
        }

        left_bin_overlap = (loc_fine[fine_start+1,2] - bin_start) / (loc_fine[fine_start+1,2] - loc_fine[fine_start,2])
        bins_left = (fine_mid - fine_start - 1) + left_bin_overlap
        right_bin_overlap = (bin_end - loc_fine[fine_end-1,2]) / (loc_fine[fine_end,2] - loc_fine[fine_end-1,2])
        bins_right = (fine_end - fine_mid - 1) + right_bin_overlap
        #wgt_left = bins_left / (bins_left+bins_right)
        #wgt_right = bins_right / (bins_left+bins_right)

        curr_wgt = bin_wgts[bin_id]
        #bin_wgts[bin_id] = curr_wgt/wgt_left
        #bin_wgts = append(bin_wgts, curr_wgt/wgt_right, after=bin_id)
        bin_wgts[bin_id] = fine_mult / bins_left
        bin_wgts = append(bin_wgts, fine_mult / bins_right, after=bin_id)
        #loc = append(loc, loc_fine[fine_mid,], after=bin_id-1)
        loc = rbind(loc[1:(bin_id-1),], loc_fine[fine_mid,], loc[bin_id:length(loc[,1]),])

        bin_raw = raw[bin_id]

        left_raw = sum(raw_fine[(fine_start+2):fine_mid,k]) + raw_fine[(fine_start+1),k] * left_bin_overlap
        right_raw = sum(raw_fine[(fine_mid+1):(fine_end-1),k]) + raw_fine[fine_end,k] * right_bin_overlap
        raw[bin_id] = left_raw

        raw = append(raw, right_raw, after=bin_id)

        #GC = rbind(GC[1:(bin_id+i-1),], GC[(bin_id+i-1),], GC[(bin_id+i-1):length(GC[,1]),])
        GC = append(GC, GC[(bin_id)], after=bin_id)
        F = append(F, F[(bin_id)], after=bin_id)
        for (chrom_id in 1:length(bounds[,1]))
          if (bounds[chrom_id,2] > bin_id)
            bounds[chrom_id,2] = bounds[chrom_id,2]+1
      }
      #for (i in 1:length(loc[,1]))
      #  cat(paste(loc[i,1], "\t", loc[i,2], "\n", sep=""))
      pos = cbind(c(1,bounds[,2]), c(bounds[,2], l))

      for (raw_i in 1:length(raw))
        raw[raw_i] = raw[raw_i] * bin_wgts[raw_i]

      normal = lowess.gc( GC, (raw+1)/mean(raw+1) )
    }
  }

  # ----------------------------------------------------------------------------
  # -- Segment data
  # ----------------------------------------------------------------------------

  # Compute log ratio between kth sample and reference
  if (stat == 0)
    lr = log2(normal)
  else
    lr = log2((normal)/(F))

  t2 <- Sys.time()
  total_other <<- total_other + (t2-t1)

  # Determine breakpoints and extract chrom/locations
  CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
  frag         = segs$output[,2:3]

  t3 <- Sys.time()
  total_seg <<- total_seg + (t3-t2)

  # Map breakpoints to kth sample
  len = dim(frag)[1]
  bps = array(0, len)
  for (j in 1:len)
    bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
  bps = sort(bps)
  bps[(len=len+1)] = l

  # Track global breakpoint locations
  seg_breaks = matrix(0,length(normal),1)
  seg_breaks[bps] = 1

  # Find outlier bins to subdivide;
  # Modify bins to contain median read count/bin within each segment
  fixed = matrix(0,length(normal),1)
  last_cn = median(normal[1:bps[2]])
  last_mean = mean(normal[1:bps[2]])
  last_sd = sd(normal[1:bps[2]])
  last_raw_mean = mean(raw[1:bps[2]])
  last_raw_sd = sd(raw[1:bps[2]])
  fixed[1:bps[2]] = last_cn

  outliers = c()
  breaks = c()
  j = 1
  chrom_i = 1

  zthresh = 1

  for(i in 2:(len-1)) {
    curr_cn = median(normal[bps[i]:(bps[i+1]-1)])
    curr_mean = mean(normal[bps[i]:(bps[i+1]-1)])
    curr_sd = sd(normal[bps[i]:(bps[i+1]-1)])
    curr_raw_mean = mean(raw[bps[i]:(bps[i+1]-1)])
    curr_raw_sd = sd(raw[bps[i]:(bps[i+1]-1)])
    fixed[bps[i]:(bps[i+1]-1)] = curr_cn

    cn_change = curr_cn - last_cn

    if (loc[bps[i]-1,1] != loc[bps[i],1]) {
      next
    }

    num_b = length(bounds[,1])
    while (chrom_i <= num_b && bounds[chrom_i,2] < bps[i])
      chrom_i = chrom_i + 1

    if ((chrom_i <= num_b) && (bounds[chrom_i,2] != bps[i]) && (abs(cn_change) >= cn_change_threshold)) {
      z1 = (normal[bps[i]-1] - last_mean) / last_sd
      z2 = (normal[bps[i]] - curr_mean) / curr_sd
      found_outlier = F

      if (curr_cn > last_cn) {
        if ((z1 > zthresh) && (z2 > 0 || z1 > -z2)) {
          if (split_mode == 3) {
            outliers[j] = bps[i]
          } else {
            outliers[j] = bps[i]-1
          }
          found_outlier = T
        }
        else if (z2 < -zthresh) {
          outliers[j] = bps[i]
          found_outlier = T
        }
      }
      else {
        if ((z1 < -zthresh) && (z2 < 0 || -z1 > z2)) {
          if (split_mode == 3) {
            outliers[j] = bps[i]
          } else {
            outliers[j] = bps[i]-1
          }
          found_outlier = T
        }
        else if (z2 > zthresh) {
          outliers[j] = bps[i]
          found_outlier = T
        }
      }

      if (split_mode == 3) {
        found_outlier = T
        outliers[j] = bps[i]
      }
      if (found_outlier) {
        if (split_mode == 0) {
          # Split in half
          breaks[j] = 0.5
        }
        else if (split_mode == 1) {
          # Simple interpolation
          #breaks[j] = (normal[outliers[j]] - last_mean) / (curr_mean - last_mean)
          breaks[j] = (curr_mean - normal[outliers[j]]) / (curr_mean - last_mean)
          if (breaks[j] < 0) {
            #cat(paste(j, ", ", breaks[j], " - ", last_mean, ", ", curr_mean, ", ", normal[bps[i]-1], ", ", normal[bps[i]], ", ", normal[outliers[j]], "\n", sep=""))
            breaks[j] = 0
            #quit()
          } else if (breaks[j] > 1) {
            #cat(paste(j, ", ", breaks[j], " - ", last_mean, ", ", curr_mean, ", ", normal[bps[i]-1], ", ", normal[bps[i]], ", ", normal[outliers[j]], "\n", sep=""))
            breaks[j] = 1
            #quit()
          }
        }
        else if (split_mode == 2) {
          # z-score-based interpolation
          zscore1 = abs((normal[outliers[j]] - last_mean)) / last_sd
          zscore2 = abs((normal[outliers[j]] - curr_mean)) / curr_sd
          breaks[j] = zscore2 / (zscore1 + zscore2)
        }
        else if (improveBounds && split_mode == 3) {
          bin_id = outliers[j]
          bin_start = loc[bin_id-2,2]
          bin_end = loc[bin_id,2]

          fine_start = 1
          found_chrom = FALSE
          for (chrom_id in 1:length(bounds_fine[,1])) {
            if (bounds_fine[chrom_id,1] == loc[bin_id,1]) {
              found_chrom = TRUE
              if (chrom_id == 1)
                fine_start = 1
              else
                fine_start = bounds_fine[chrom_id-1,2]
              break
            }
          }
          if (!found_chrom) {
            cat(paste("Couldn't find chromosome ", loc[bin_id,1], "\n", sep=""))
            quit()
          }
          while (loc_fine[fine_start,1] != loc[bin_id,1] || loc_fine[fine_start,2] <= bin_start) # Find the right bin
            fine_start = fine_start+1

          fine_end = fine_start
          while (loc_fine[fine_end+1,2] < bin_end)
            fine_end = fine_end+1

          left_bin_overlap = (loc_fine[fine_start,2] - bin_start) / (loc_fine[fine_start,2] - loc_fine[fine_start-1,2])
          right_bin_overlap = (bin_end - loc_fine[fine_end,2]) / (loc_fine[fine_end+1,2] - loc_fine[fine_end,2])

          total_bins = fine_end - fine_start
          total_bins_olap = total_bins + left_bin_overlap + right_bin_overlap
          left_bins = 0
          right_bins = total_bins
          left_count = 0
          right_count = sum(raw_fine[(fine_start+1):fine_end,k])

          #cat(paste("Last: ", last_raw_mean, " (", last_raw_sd, ")", "\n", sep=""))
          #cat(paste("Next: ", curr_raw_mean, " (", curr_raw_sd, ")", "\n", sep=""))
          #cat(paste("Raw: ", paste(raw_fine[fine_start:fine_end,k]*total_bins_olap/2, collapse=","), "\n", sep=""))

          #last_ln_sd = log(last_raw_sd)
          #curr_ln_sd = log(curr_raw_sd)
          #prob_left = c()
          #prob_right = c()
          #for (curr_bin_id in 1:total_bins) {
            #z_left = (raw_fine[fine_start+curr_bin_id,k]*total_bins_olap/2 - last_raw_mean) / last_raw_sd
            #z_right = (raw_fine[fine_start+curr_bin_id,k]*total_bins_olap/2 - curr_raw_mean) / curr_raw_sd
            #prob_left[curr_bin_id] = (z_left ** 2)/2 + last_ln_sd
            #prob_right[curr_bin_id] = (z_right ** 2)/2 + curr_ln_sd
          #  sq_left = (raw_fine[fine_start+curr_bin_id,k]*total_bins_olap/2 - last_raw_mean) ** 2
          #  sq_right = (raw_fine[fine_start+curr_bin_id,k]*total_bins_olap/2 - curr_raw_mean) ** 2
          #  prob_left[curr_bin_id] = sq_left
          #  prob_right[curr_bin_id] = sq_right
          #}

          #cat(paste("Left p's: ", paste(prob_left, collapse=","), "\n", sep=""))
          #cat(paste("Right p's: ", paste(prob_right, collapse=","), "\n", sep=""))

          #cat(paste(sum(prob_left[1:100]), " - ", sum(prob_left[101:total_bins]), "\n", sep=""))
          #cat(paste(sum(prob_right[1:100]), " - ", sum(prob_right[101:total_bins]), "\n", sep=""))
          #cat(paste(paste(prob_left, collapse=","), "\n", sep=""))
          #cat(paste(paste(prob_right, collapse=","), "\n", sep=""))

          #curr_score = sum(prob_right)
          #min_score = curr_score
          #min_score = abs(right_count*total_bins_olap/(2*right_bins) - curr_raw_mean) / curr_raw_sd
          #if (curr_raw_mean > last_raw_mean)
          #  min_score = -min_score
          min_score = (right_count*fine_mult/right_bins - curr_raw_mean) ** 2
          peak_left_bins = left_bins

          margin = 0
          for (left_bins in 1:(total_bins-margin)) {
            right_bins = total_bins - left_bins
            left_count = left_count + raw_fine[fine_start+left_bins,k]
            right_count = right_count - raw_fine[fine_start+left_bins,k]
            #zscore1 = abs(left_count*total_bins_olap/(2*left_bins) - last_raw_mean) / last_raw_sd
            #if (last_raw_mean < curr_raw_mean)
            #  zscore1 = -zscore1
            #if (zscore1 < 0)
            #  zscore1 = 0
            zscore1 = (left_count*fine_mult/left_bins - last_raw_mean) ** 2
            curr_score = zscore1
            if (right_bins > 0) {
                #zscore2 = abs(right_count*total_bins_olap/(2*right_bins) - curr_raw_mean) / curr_raw_sd
                #if (curr_raw_mean < last_raw_mean)
                #  zscore2 = -zscore2
                #if (zscore2 < 0)
                #  zscore2 = 0
                zscore2 = (right_count*fine_mult/right_bins - curr_raw_mean) ** 2
                curr_score = (zscore1 + zscore2)/2
            #    if (outliers[j] == 2904) {
            #      cat(paste(left_bins, " left (count = ", left_count, ", z = ", zscore1, ")\t", right_bins, " right (count = ", right_count, ", z = ", zscore2, ")\tScore = ", curr_score, "\n", sep=""))
            #    }
            } #else {
            #    if (outliers[j] == 2904) {
            #      cat(paste(left_bins, " left (count = ", left_count, ", z = ", zscore1, ")\t", right_bins, " right (count = 0, z = 0)\t\tScore = ", curr_score, "\n", sep=""))
            #    }
            #}

            #curr_score = curr_score - prob_right[left_bins] + prob_left[left_bins]

            if (margin == left_bins || curr_score < min_score) {
                min_score = curr_score
                peak_left_bins = left_bins
            }
          }
          new_bound = loc_fine[fine_start+peak_left_bins,2]
          #if (new_bound < loc[bin_id-1,2]) {
          loc[bin_id-1,2] = new_bound
          #} else {
          #  loc[bin_id,2] = new_bound
          #}
        }
        j = j+1
      }
    }

    last_cn = curr_cn
    last_mean = curr_mean
    last_sd = curr_sd
    last_raw_mean = curr_raw_mean
    last_raw_sd = curr_raw_sd
  }
  normal = normal/mean(fixed)
  fixed = fixed/mean(fixed)

  # ----------------------------------------------------------------------------
  # -- Determine Copy Number (SoS Method)
  # ----------------------------------------------------------------------------

  # Determine Copy Number     
  CNgrid           = seq(1.5, maxPloidy, by=0.05)
  outerRaw         = fixed %o% CNgrid
  outerRound       = round(outerRaw)
  outerDiff        = (outerRaw - outerRound) ^ 2
  outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult[,k]       = CNgrid[order(outerColsums[,k])]
  CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

  if (f == 0 | length(which(lab[k]==ploidy[,1]))==0 ) {
    CN = CNmult[1,k]
  } else if (f == 1) {
    CN = ploidy[which(lab[k]==ploidy[,1]),2]
  } else {
    estimate = ploidy[which(lab[k]==ploidy[,1]),2]
    CN = CNmult[which(abs(CNmult[,k] - estimate)<.4),k][1]
  }
  final = round(fixed*CN)

  # Output results of CN calculations to file
  out=paste(lab[k], CN, paste(CNmult[,k], collapse= ","), paste(CNerror[,k], collapse= ","), sep="\t")
  cat(out, "\n")

  #if (length(bin_wgts) > 0) {
  #  wgts_col = matrix(0,length(bin_wgts),1)
  #  for (index in 1:length(bin_wgts))
  #    wgts_col[index,1] = bin_wgts[index]

  #  #wgts_col[1,] = bin_wgts
  cat(paste(length(raw), "\n", sep=""))
  write.table(cbind(loc,final,raw), file=paste(lab[k], "_CN", plot_sub, ".tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc),lab[k], "Raw"), sep="\t", quote=FALSE)
  #}
  #else {
  #write.table(cbind(loc,final), file=paste(lab[k], "_CN", plot_sub, ".tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc),lab[k]), sep="\t", quote=FALSE)
  #}

  t4 <- Sys.time()
  total_other <<- total_other + (t4-t3)

  if (missing(plot_sub)) {
    if (length(raw) == length(bin_wgts)) {
      for (raw_i in 1:length(raw))
        raw[raw_i] = raw[raw_i] / bin_wgts[raw_i]
    }

    v <- list("out" = outliers, "breaks" = breaks, "wgt" = bin_wgts, "loc" = loc, "gc" = GC, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
    return(v)
  }

  # ----------------------------------------------------------------------------
  # -- Generate Plots & Figures
  # ----------------------------------------------------------------------------

  tPlotA <- proc.time()

  # Plot Distribution of Read Coverage
  jpeg(filename=paste(lab[k], "_dist", plot_sub, ".jpeg", sep=""), width=3000, height=750)
  
  top=round(quantile(raw, c(.995))[[1]])
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  main=data.frame(x=which(raw<top), y=raw[which(raw<top)])
  dist_outliers=data.frame(x=which(raw>top), y=array(top*.99, length(which(raw>top))))
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))

  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) + 
    geom_point(data=main, aes(x=x, y=y), size=3) +
    geom_point(data=dist_outliers, aes(x=x, y=y), shape=5, size=4) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    labs(title=paste("Genome Wide Read Distribution for Sample \"", lab[k], "\"", sep=""), x="Chromosome", y="Read Count", size=16) +
    theme(plot.title=element_text(size=36, vjust=1.5)) +
    theme(axis.title.x=element_text(size=40, vjust=-.1), axis.title.y=element_text(size=40, vjust=-.06)) +
    theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'gray90')) +
    theme(plot.margin=unit(c(.5,1,.5,1.5),"cm")) +
    theme(panel.grid.major.x = element_blank()) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    geom_vline(xintercept = c(1, l), size=.25) +
    geom_hline(yintercept = c(-top*.1, top), size=.25)

    grid.arrange(plot1, ncol=1)
  dev.off()

  tPlotB <- proc.time()
  cat(paste("    Dist: ", (tPlotB-tPlotA)['elapsed'], "\n", sep=""))

  #Plot histogram of bin counts
  jpeg(filename=paste(lab[k], "_counts", plot_sub, ".jpeg", sep=""), width=2500, height=1500)
    par(mar = c(7.0, 7.0, 7.0, 3.0))

    temp=sort(raw)[round(l*.01) : (l-round(l*.01))] 
    reads = hist(temp, breaks=100, plot=FALSE)
    # plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", xaxt="n", cex.main=3, cex.axis=2, cex.lab=2)
    plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", cex.main=3, cex.axis=2, cex.lab=2)
    # axis(side=1, at=seq(min(temp), round(diff(range(temp))/20)*22, round(diff(range(temp))/20)), cex.axis=2)
    tu = par('usr')
    par(xpd=FALSE)
    clip(tu[1], as.integer(mean(temp)-(diff(reads$mids)/2)), tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(mean(temp)+(diff(reads$mids)/2), tu[2], tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(tu[1], mean(temp) - sd(temp), tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(mean(temp) + sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(tu[1], mean(temp) - 2*sd(temp), tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE)
    clip(mean(temp) + 2*sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE) 
    legend("topright", inset=.05, legend=c("mean", "< 1σ", "> 1σ", "> 2σ"), fill=c("black", "gray50", "gray75", "gray90"), cex=2.5)
  dev.off()

  tPlotC <- proc.time()
  cat(paste("    Counts: ", (tPlotC-tPlotB)['elapsed'], "\n", sep=""))

  #Plot lorenz curves
  jpeg(filename=paste(lab[k], "_lorenz", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  nReads=sum(raw)
  uniq=unique(sort(raw))
  
  lorenz=matrix(0, nrow=length(uniq), ncol=2)
  a=c(length(which(raw==0)), tabulate(raw, nbins=max(raw)))
  b=a*(0:(length(a)-1))
  for (i in 2:length(uniq)) {
    lorenz[i,1]=sum(a[1:uniq[i]])/l
    lorenz[i,2]=sum(b[2:uniq[i]])/nReads
  }

  # smooth.spline needs >= 4 points...
  fit = data.frame(x=lorenz[,1], y=lorenz[,2])
  if(nrow(lorenz) >= 4)
  {
    spline = try(smooth.spline(lorenz))
    if(class(spline) != "try-error")
      fit = data.frame(x=spline$x, y=spline$y)
  }

  perf=data.frame(x=c(0,1), y=c(0,1))

  plot1 = try(ggplot() +
    geom_line(data=perf, aes(x=x, y=y, color="Perfect Uniformity"), size=3) +
    geom_line(data=fit, aes(x=x, y=y, color="Sample Uniformity"), size=3) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    labs(title=paste("Lorenz Curve of Coverage Uniformity for Sample ", lab[k], sep=""), x="Cumulative Fraction of Genome", y="Cumulative Fraction of Total Reads") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) + theme(legend.title=element_blank(), legend.text=element_text(size=40)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.15, .85)) +
    scale_color_manual(name='', values=c('Perfect Uniformity'="black", 'Sample Uniformity'=colors[cp,1])))

    grid.arrange(plot1, ncol=1)
  dev.off()

  tPlotD <- proc.time()
  cat(paste("    Lorenz: ", (tPlotD-tPlotC)['elapsed'], "\n", sep=""))

  ##Plot GC correction
  jpeg(filename=paste(lab[k], "_GC", outSuffix, ".jpeg", sep=""), width=2500, height=1250)

  low = lowess(GC, log(normal2[,k]), f=0.05)
  app = approx(low$x, low$y, GC)
  cor = exp(log(normal2[,k]) - app$y)
  
  uncorrected = data.frame(x=GC, y=log(normal2[,k]))
  corrected = data.frame(x=GC, y=log(cor))
  fit = data.frame(x=app$x, y=app$y)

  try(plot1 <- ggplot() +
    geom_point(data=uncorrected, aes(x=x, y=y), size=3) +
    geom_line(data=fit, aes(x=x, y=y, color="Lowess Fit"), size=3) +
    scale_x_continuous(limits=c(min(.3, min(GC)), max(.6, max(GC))), breaks=seq(.3,.6,.05)) +
    labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Uncorrected)", sep=""), x="GC content", y="Normalized Read Counts (log scale)") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    #theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
    #theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    #theme(legend.position=c(.85, .9)) +
    scale_color_manual(name='', values=colors[cp,1]))

  try(plot2 <- ggplot() +
    geom_point(data=corrected, aes(x=x, y=y), size=3) +
    scale_x_continuous(limits=c(min(.3, min(GC)), max(.6, max(GC))), breaks=seq(.3,.6,.05)) +
    labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Corrected)", sep=""), x="GC content", y="") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')))

    try(grid.arrange(plot1, ncol=1))
    #try(grid.arrange(plot1, plot2, ncol=2))
  dev.off()

  #Plot Scaled/Normalized Bin Count Histogram
  jpeg(filename=paste(lab[k], "_hist", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  clouds=data.frame(x=normal*CN)

  plot1 = ggplot() +
    geom_histogram(data=clouds, aes(x=x), binwidth=.05, color="black", fill="gray60") +
    geom_vline(xintercept=seq(0,10,1), size=1, linetype="dashed", color=colors[cp,1]) +
    scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1)) +
    labs(title=paste("Frequency of Bin Counts for Sample \"", lab[k], "\"\nNormalized and Scaled by Predicted CN (", CNmult[1,k], ")", sep=""), x="Copy Number", y="Frequency") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black'))

    grid.arrange(plot1, ncol=1)
  dev.off()

  tPlotE <- proc.time()
  cat(paste("    Hist: ", (tPlotE-tPlotD)['elapsed'], "\n", sep=""))

  #Plot sum of squares error for each potential copy number
  jpeg(filename=paste(lab[k], "_SoS", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  top = max(outerColsums[,k])
  dat = data.frame(x=CNgrid, y=outerColsums[,k])
  lim = cbind(c(seq(0,5000,500), 1000000), c(50, 100, 100, 200, 250, 400, 500, 500, 600, 600, 750, 1000))
  step = lim[which(top<lim[,1])[1],]
  minSoS = data.frame(x=CNmult[1,k], y=CNerror[1,k])
  bestSoS = data.frame(x=CN, y=outerColsums[which(CNgrid==CN),k])

  plot1 = ggplot() +
    geom_line(data=dat, aes(x=x, y=y), size=3) +
    geom_point(data=dat, aes(x=x, y=y), shape=21, fill="black", size=5) +
    geom_point(data=minSoS, aes(x=x, y=y*1.02, color="Minimum SoS Error"), shape=18, size=15) +
    geom_point(data=bestSoS, aes(x=x, y=y*.98, color="Chosen Ploidy"), shape=18, size=15) +
    scale_x_continuous(limits=c(1.5, 6), breaks=seq(1.5, 6, .5)) +
    scale_y_continuous(limits=c(.5*min(outerColsums[,k]), top), breaks=seq(0, step[1], step[2])) +
    labs(title="Sum of Squares Error Across Potential Copy Number States", x="Copy Number Multiplier", y="Sum of Squares Error") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.85, .9)) +
    scale_color_manual(name='', values=c('Minimum SoS Error'=colors[cp,1], 'Chosen Ploidy'=colors[cp,2]))

    grid.arrange(plot1, ncol=1)
  dev.off()

  tPlotF <- proc.time()
  cat(paste("    SoS: ", (tPlotF-tPlotE)['elapsed'], "\n", sep=""))

  #Plot colored CN profile

  top=8
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  clouds=data.frame(x=1:length(normal), y=normal*CN)
  amp=data.frame(x=which(final>2), y=final[which(final>2)])
  del=data.frame(x=which(final<2), y=final[which(final<2)])
  flat=data.frame(x=which(final==2), y=final[which(final==2)])
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))
  outlier_pts=data.frame(x=outliers, y=normal[outliers]*CN)

  jpeg(filename=paste(lab[k], "_CN", plot_sub, ".jpeg", sep=""), width=3000, height=750)
  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
    geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
    geom_point(data=flat, aes(x=x, y=y), size=4) +
    geom_point(data=amp, aes(x=x, y=y), size=4, color='black') +
    geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
    geom_point(data=outlier_pts, aes(x=x, y=y), size=5, color=colors[cp,1]) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n Predicted Ploidy = ", CN, sep=""), x="Chromosome", y="Copy Number", size=16) +
    theme(plot.title=element_text(size=40, vjust=1.5)) +
    theme(axis.title.x=element_text(size=40, vjust=-.05), axis.title.y=element_text(size=40, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'gray90')) +
    theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
    theme(panel.grid.major.x = element_blank()) +
    geom_vline(xintercept = c(1, l), size=.5) +
    geom_hline(yintercept = c(-top*.1, top), size=.5) +
    xlim(pos[7]+60,pos[8]-1)

    grid.arrange(plot1, ncol=1)
  dev.off()

  tPlotG <- proc.time()
  cat(paste("    CN: ", (tPlotG-tPlotF)['elapsed'], "\n", sep=""))

  t5 <- Sys.time()
  total_plots <<- total_plots + (t5-t4)

  if (length(raw) == length(bin_wgts)) {
    for (raw_i in 1:length(raw))
      raw[raw_i] = raw[raw_i] / bin_wgts[raw_i]
  }
  v <- list("out" = outliers, "breaks" = breaks, "wgt" = bin_wgts, "loc" = loc, "gc" = GC, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
  return(v)
}

processSampleNoNewBins <- function(k, loc, raw, GC, bounds, outliers, breaks, iter=1, plot_sub=NULL) {
  statusFile = file( paste(user_dir, "/", status, sep="") )
  writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>", lab[k], "</processingfile>", sep=""), paste("<percentdone>", (k*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
  close(statusFile)

  t1 <- Sys.time()

  if (missing(plot_sub))
    cat(paste("Processing sample ", lab[k], "\n", sep=""))
  else
    cat(paste("Processing sample ", lab[k], plot_sub, "\n", sep=""))

  # ----------------------------------------------------------------------------
  # -- Split outlier bins
  # ----------------------------------------------------------------------------
  cat("Splitting outlier bins\n")
  #cat(paste(length(loc[,1]), " bins\n", sep=''))
  #normal = normal[,k]
  if (!missing(outliers)) {
    num_o = length(outliers)
    if (num_o > 0) {
      cat(paste(length(outliers), " outliers: ", paste(outliers, collapse=','), "\n", sep=''))
      cat(paste("Breaks: ", paste(breaks, collapse=','), "\n", sep=''))
      bin_wids = c()
      for (i in 1:length(raw))
        bin_wids[i] = 1

      for (i in 1:num_o) {
        #otimeA <- proc.time()

        bin_id = outliers[i]
        #cat(paste(i, " / ", num_o, ", ", bin_id, ",", breaks[i], "\n", sep=''))
        #cat(paste(paste(loc[bin_id,], collapse=','), "\n", sep=''))
        if (loc[bin_id-1,1] != loc[bin_id,1])
          next 

        bin_start = loc[bin_id-1,2]
        bin_end = loc[bin_id,2]
        bin_mid = (bin_end - bin_start) * breaks[i] + bin_start

        #otimeB <- proc.time()

        fine_start = 1
        found_chrom = FALSE
        for (chrom_id in 1:length(bounds_fine[,1])) {
          if (bounds_fine[chrom_id,1] == loc[bin_id,1]) {
            found_chrom = TRUE
            if (chrom_id == 1)
              fine_start = 1
            else
              fine_start = bounds_fine[chrom_id-1,2]
            break
          }
        }
        if (!found_chrom) {
          cat(paste("Couldn't find chromosome ", loc[bin_id,1], "\n", sep=""))
          quit()
        }

        #cat(paste("fine_start = ", fine_start, "\n", sep=""))
        #cat(paste("bin_start = ", loc[bin_id-1,1], " - ", bin_start, "\n", sep=""))
        while (loc_fine[fine_start+1,2] <= bin_start) # Find the right bin
          fine_start = fine_start+1

        #otimeC <- proc.time()

        fine_mid = fine_start
        while (loc_fine[fine_mid,2] < bin_mid)
          fine_mid = fine_mid+1
        if ((bin_mid - loc_fine[fine_mid-1,2]) < (loc_fine[fine_mid,2] - bin_mid))
          fine_mid = fine_mid - 1

        #otimeD <- proc.time()

        fine_end = fine_start
        while (loc_fine[fine_end,2] < bin_end)
          fine_end = fine_end+1

        #otimeE <- proc.time()

        left_bin_overlap = (loc_fine[fine_start+1,2] - bin_start) / (loc_fine[fine_start+1,2] - loc_fine[fine_start,2])
        bins_left = (fine_mid - fine_start - 1) + left_bin_overlap
        right_bin_overlap = (bin_end - loc_fine[fine_end-1,2]) / (loc_fine[fine_end,2] - loc_fine[fine_end-1,2])
        bins_right = (fine_end - fine_mid - 1) + right_bin_overlap
        wgt_left = bins_left / (bins_left+bins_right)
        wgt_right = bins_right / (bins_left+bins_right)
        left_raw = sum(raw_fine[(fine_start+2):fine_mid,k]) + raw_fine[(fine_start+1),k] * left_bin_overlap
        right_raw = sum(raw_fine[(fine_mid+1):(fine_end-1),k]) + raw_fine[fine_end,k] * right_bin_overlap

        if (bin_wids[bin_id-1] != 1) {
          cat(paste("Last bin was also split!", "\n", sep=""))
          quit()
        }

        if (breaks[i] < 0.5) {
          bin_wids[bin_id-1] = bin_wids[bin_id-1] + wgt_left
          bin_wids[bin_id] = wgt_right
          loc[bin_id-1,] = loc_fine[fine_mid,]
          raw[bin_id-1] = (raw[bin_id-1] + left_raw) / (1 + bins_left/fine_mult)
          raw[bin_id] = right_raw * fine_mult / bins_right
          #raw[bin_id-1] = (raw[bin_id-1] + left_raw) / bin_wids[bin_id-1]
          #raw[bin_id] = right_raw / bin_wids[bin_id]
        }
        else {
          bin_wids[bin_id] = wgt_left
          bin_wids[bin_id+1] = bin_wids[bin_id+1] + wgt_right
          loc[bin_id,] = loc_fine[fine_mid,]
          raw[bin_id] = left_raw * fine_mult / bins_left
          raw[bin_id+1] = (raw[bin_id+1] + right_raw) / (1 + bins_right/fine_mult)
          #raw[bin_id] = left_raw / bin_wids[bin_id]
          #raw[bin_id+1] = (raw[bin_id+1] + right_raw) / bin_wids[bin_id+1]
        }

        #otimeF <- proc.time()
        #cat(paste((otimeF-otimeA)['elapsed'], "\t:\t", (otimeB-otimeA)['elapsed'], "\t", (otimeC-otimeB)['elapsed'], "\t", (otimeD-otimeC)['elapsed'], "\t", (otimeE-otimeD)['elapsed'], "\t", (otimeF-otimeE)['elapsed'], "\n", sep=""))
      }
      pos = cbind(c(1,bounds[,2]), c(bounds[,2], l))

      normal[,k] = lowess.gc( GC, (raw+1)/mean(raw+1) )
    }
  }
  
  cat("Finished splitting outlier bins\n")

  #cat("Normal:\n")
  #cat(paste(paste(normal,collapse=','), "\n", sep=''))
  #cat("Segmenting\n")

  # ----------------------------------------------------------------------------
  # -- Segment data
  # ----------------------------------------------------------------------------

  # Compute log ratio between kth sample and reference
  if (stat == 0)
    lr = log2(normal[,k])
  else
    lr = log2((normal[,k])/(F))

  t2 <- Sys.time()
  total_other <<- total_other + (t2-t1)

  # Determine breakpoints and extract chrom/locations
  CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
  frag         = segs$output[,2:3]

  #cat("Finished segmenting\n")

  t3 <- Sys.time()
  total_seg <<- total_seg + (t3-t2)

  # Map breakpoints to kth sample
  len = dim(frag)[1]
  bps = array(0, len)
  for (j in 1:len)
    bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
  bps = sort(bps)
  bps[(len=len+1)] = l

  #cat(paste(length(loc,1), "\n", sep=''))
  #cat(paste(paste(bps, collapse=','), sep=''))

  # Track global breakpoint locations
  seg_breaks = matrix(0,l,1)
  seg_breaks[bps] = 1

  #cat("Finding outlier bins\n")

  cat(paste(paste(bps, collapse=","), "\n", sep=""))

  # Find outlier bins to subdivide;
  # Modify bins to contain median read count/bin within each segment
  fixed = matrix(0,l,1)
  last_cn = median(normal[1:bps[2],k])
  last_mean = mean(normal[1:bps[2],k])
  last_sd = sd(normal[1:bps[2],k])
  fixed[1:bps[2]] = last_cn

  new_outliers = c()
  new_breaks = c()
  j = 1
  chrom_i = 1

  zthresh = 1

  for(i in 2:(len-1)) {
    cat(paste("BP: ", bps[i], "\n", sep=""))
    curr_cn = median(normal[bps[i]:(bps[i+1]-1),k])
    curr_mean = mean(normal[bps[i]:(bps[i+1]-1),k])
    curr_sd = sd(normal[bps[i]:(bps[i+1]-1),k])
    fixed[bps[i]:(bps[i+1]-1)] = curr_cn

    cn_change = curr_cn - last_cn
    num_b = length(bounds[,1])
    while (chrom_i <= num_b && bounds[chrom_i,2] < bps[i])
      chrom_i = chrom_i + 1

    if (is.na(cn_change)) {
      next
    }

    #cat(paste(chrom_i, ", ", num_b, ", ", cn_change, ", ", cn_change_threshold, "\n", sep=""))

    if (!is.na(normal[bps[i]-1,k]) && !is.na(normal[bps[i],k]) && (chrom_i <= num_b) && (bounds[chrom_i,2] != bps[i]) && (abs(cn_change) >= cn_change_threshold)) {
      z1 = (normal[bps[i]-1,k] - last_mean) / last_sd
      z2 = (normal[bps[i],k] - curr_mean) / curr_sd
      cat(paste(last_mean, " (", last_sd, "), ", curr_mean, " (", curr_sd, ") - ", z1, ", ", z2, "\n", sep=""))
      found_outlier = F
      if (curr_cn > last_cn) {
        if ((z1 > zthresh) && (z2 > 0 || z1 > -z2)) {
          new_outliers[j] = bps[i]-1
          found_outlier = T
        }
        else if (z2 < -zthresh) {
          new_outliers[j] = bps[i]
          found_outlier = T
        }
      }
      else {
        if ((z1 < -zthresh) && (z2 < 0 || -z1 > z2)) {
          new_outliers[j] = bps[i]-1
          found_outlier = T
        }
        else if (z2 > zthresh) {
          new_outliers[j] = bps[i]
          found_outlier = T
        }
      }
      if (found_outlier) {
        if (split_mode == 0) {
          # Split in half
          new_bin_wid = 0.5 ** iter
          cat(paste(bin_wids[new_outliers[j]-1], ", ", bin_wids[new_outliers[j]], ", ", bin_wids[new_outliers[j]+1], "\n", sep=""))
          if (bin_wids[new_outliers[j]] > 1) {
            if (bin_wids[new_outliers[j]+1] < 1) {
              new_breaks[j] = bin_wids[new_outliers[j]] - 1 - new_bin_wid
              new_outliers[j] = new_outliers[j]+1
            }
            else {
              new_breaks[j] = bin_wids[new_outliers[j]-1] + new_bin_wid
              new_outliers[j] = new_outliers[j]-1
            }
          }
          else if (bin_wids[new_outliers[j]] < 1) {
            if (bin_wids[new_outliers[j]+1] > 1) {
              new_breaks[j] = bin_wids[new_outliers[j]] + new_bin_wid
            }
            else {
              new_breaks[j] = 1 - bin_wids[new_outliers[j]] - new_bin_wid
            }
          }
          else {
            new_breaks[j] = 0.5
          }
          #quit()
        }
        else {
          if (split_mode == 1) {
            # Simple interpolation
            new_breaks[j] = (curr_mean - normal[new_outliers[j],k]) / (curr_mean - last_mean)
          }
          else if (split_mode == 2) {
            # z-score-based interpolation
            zscore1 = abs((normal[new_outliers[j],k] - last_mean)) / last_sd
            zscore2 = abs((normal[new_outliers[j],k] - curr_mean)) / curr_sd
            new_breaks[j] = zscore2 / (zscore1 + zscore2)
          }

          # Scale offset in bin by bin width
          scaled_break = new_breaks[j] * bin_wids[new_outliers[j]]
          # Shift beginning of bin to original bin boundaries
          prev_bin_wid = bin_wids[new_outliers[j]-1]
          scaled_break = scaled_break + (prev_bin_wid + bin_wids[new_outliers[j]-2] - 2)
        
          cat(paste(new_outliers[j], ", ", scaled_break, "\n", sep=""))

          if (scaled_break < 0) {
            cat("Found one < 0\n")
            new_outliers[j] = new_outliers[j]-1
            new_breaks[j] = 1 + scaled_break
          }
          else if (scaled_break > 1) {
            cat("Found one > 1\n")
            new_outliers[j] = new_outliers[j]+1
            new_breaks[j] = scaled_break - 1
          }
          else {
            new_breaks[j] = scaled_break
          }
        }

        j = j+1
      }
    #for (o_id in 1:num_o) {
    #  if outliers[o_id] == bps[i]
        
    }

    last_cn = curr_cn
    last_mean = curr_mean
    last_sd = curr_sd
  }

  cat(paste("Old outliers: ", paste(outliers, collapse=","), "\n", sep=""))
  cat(paste("Old Breaks:   ", paste(breaks, collapse=","), "\n", sep=""))
  cat(paste("New outliers: ", paste(new_outliers, collapse=","), "\n", sep=""))
  cat(paste("New Breaks:   ", paste(new_breaks, collapse=","), "\n", sep=""))

  num_prev_o = length(outliers)
  num_curr_o = length(new_outliers)
  if (num_curr_o > 0) {
    curr_id = 1
    for (prev_id in 1:num_prev_o) {
      #cat(paste("  ", prev_id, "\n", sep=""))
      while (curr_id <= num_curr_o && (outliers[prev_id] > (new_outliers[curr_id]+2)))
        curr_id = curr_id+1

      if ((curr_id > num_curr_o) || (outliers[prev_id] < (new_outliers[curr_id]-2))) {
        if (curr_id == 1) {
          new_outliers = c(outliers[prev_id], new_outliers[curr_id:num_curr_o])
          new_breaks = c(breaks[prev_id], new_breaks[curr_id:num_curr_o])
        }
        else if (curr_id > num_curr_o) {
          new_outliers = c(new_outliers[1:(curr_id-1)], outliers[prev_id])
          new_breaks = c(new_breaks[1:(curr_id-1)], breaks[prev_id])
        }
        else {
          new_outliers = c(new_outliers[1:(curr_id-1)], outliers[prev_id], new_outliers[curr_id:num_curr_o])
          new_breaks = c(new_breaks[1:(curr_id-1)], breaks[prev_id], new_breaks[curr_id:num_curr_o])
        }
        num_curr_o = num_curr_o+1
      }
      cat(paste("New outliers: ", paste(new_outliers, collapse=","), "\n", sep=""))
    }
    cat(paste("Combined:     ", paste(new_outliers, collapse=","), "\n", sep=""))
    cat(paste("Breaks:       ", paste(new_breaks, collapse=","), "\n", sep=""))
    outliers = new_outliers
    breaks = new_breaks
  }
  #cat(paste("\nOutliers: ", paste(outliers, collapse=","), "\n", sep=''))
  normal[,k] = normal[,k]/mean(fixed)
  fixed = fixed/mean(fixed)

  # ----------------------------------------------------------------------------
  # -- Determine Copy Number (SoS Method)
  # ----------------------------------------------------------------------------

  # Determine Copy Number     
  CNgrid           = seq(1.5, maxPloidy, by=0.05)
  outerRaw         = fixed %o% CNgrid
  outerRound       = round(outerRaw)
  outerDiff        = (outerRaw - outerRound) ^ 2
  outerColsums[,k] = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult[,k]       = CNgrid[order(outerColsums[,k])]
  CNerror[,k]      = round(sort(outerColsums[,k]), digits=2)

  if (f == 0 | length(which(lab[k]==ploidy[,1]))==0 ) {
    CN = CNmult[1,k]
  } else if (f == 1) {
    CN = ploidy[which(lab[k]==ploidy[,1]),2]
  } else {
    estimate = ploidy[which(lab[k]==ploidy[,1]),2]
    CN = CNmult[which(abs(CNmult[,k] - estimate)<.4),k][1]
  }
  final = round(fixed*CN)

  # Output results of CN calculations to file
  out=paste(lab[k], CN, paste(CNmult[,k], collapse= ","), paste(CNerror[,k], collapse= ","), sep="\t")
  cat(out, "\n")

  write.table(cbind(loc,final), file=paste(lab[k], "_CN", plot_sub, ".tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc),lab[k]), sep="\t", quote=FALSE)

  t4 <- Sys.time()
  total_other <<- total_other + (t4-t3)

  if (missing(plot_sub)) {
    v <- list("out" = outliers, "breaks" = breaks, "loc" = loc, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
    return(v)
    #return(outliers)
  }

  # ----------------------------------------------------------------------------
  # -- Generate Plots & Figures
  # ----------------------------------------------------------------------------

  #Plot colored CN profile
  top=8
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  clouds=data.frame(x=1:length(normal[,k]), y=normal[,k]*CN)
  amp=data.frame(x=which(final>2), y=final[which(final>2)])
  del=data.frame(x=which(final<2), y=final[which(final<2)])
  flat=data.frame(x=which(final==2), y=final[which(final==2)])
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))
  outlier_pts=data.frame(x=outliers, y=normal[outliers,k]*CN)

  plot_chroms <- FALSE
  if (plot_chroms) {
    for (chr_id in 1:(length(pos[,1])-1)) {
      cat(paste(pos[chr_id], ",", pos[chr_id+1], "\n", sep=''))
      jpeg(filename=paste(lab[k], "_CN", plot_sub, ".chr", chr_id, ".jpeg", sep=""), width=3000, height=750)
      plot1 = ggplot() +
        geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
        geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
        geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
        geom_point(data=flat, aes(x=x, y=y), size=4) +
        geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
        geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
        geom_point(data=outlier_pts, aes(x=x, y=y), size=4, color='blue4') +
        geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
        #scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
        scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
        labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n Predicted Ploidy = ", CN, sep=""), x="Chromosome", y="Copy Number", size=16) +
        theme(plot.title=element_text(size=40, vjust=1.5)) +
        theme(axis.title.x=element_text(size=40, vjust=-.05), axis.title.y=element_text(size=40, vjust=.1)) +
        theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
        theme(panel.background = element_rect(fill = 'gray90')) +
        theme(plot.margin=unit(c(.5,1,.5,1),"cm")) +
        theme(panel.grid.major.x = element_blank()) +
        geom_vline(xintercept = c(1, l), size=.5) +
        geom_hline(yintercept = c(-top*.1, top), size=.5) +
        xlim(pos[chr_id],pos[chr_id+1])

        grid.arrange(plot1, ncol=1)
      dev.off()
    }
  }
  else {
    jpeg(filename=paste(lab[k], "_CN", plot_sub, ".jpeg", sep=""), width=3000, height=750)
    plot1 = ggplot() +
      geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
      geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
      geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
      geom_point(data=flat, aes(x=x, y=y), size=4) +
      geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
      geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
      geom_point(data=outlier_pts, aes(x=x, y=y), size=4, color='blue4') +
      geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
      scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
      scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
      labs(title=paste("Integer Copy Number Profile for Sample \"", lab[k], "\"\n Predicted Ploidy = ", CN, sep=""), x="Chromosome", y="Copy Number", size=16) +
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
  }

  # Plot Distribution of Read Coverage
  jpeg(filename=paste(lab[k], "_dist", plot_sub, ".jpeg", sep=""), width=3000, height=750)
  
  top=round(quantile(raw, c(.995))[[1]])
  rectangles1=data.frame(pos[seq(1,nrow(pos), 2),])
  rectangles2=data.frame(pos[seq(2,nrow(pos), 2),])
  main=data.frame(x=which(raw<top), y=raw[which(raw<top)])
  dist_outliers=data.frame(x=which(raw>top), y=array(top*.99, length(which(raw>top))))
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))

  plot1 = ggplot() +
    geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
    geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) + 
    geom_point(data=main, aes(x=x, y=y), size=3) +
    geom_point(data=dist_outliers, aes(x=x, y=y), shape=5, size=4) +
    geom_text(data=anno, aes(x=x, y=y, label=chrom), size=12) +
    labs(title=paste("Genome Wide Read Distribution for Sample \"", lab[k], "\"", sep=""), x="Chromosome", y="Read Count", size=16) +
    theme(plot.title=element_text(size=36, vjust=1.5)) +
    theme(axis.title.x=element_text(size=40, vjust=-.1), axis.title.y=element_text(size=40, vjust=-.06)) +
    theme(axis.text=element_text(color="black", size=40), axis.ticks=element_line(color="black"))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'gray90')) +
    theme(plot.margin=unit(c(.5,1,.5,1.5),"cm")) +
    theme(panel.grid.major.x = element_blank()) +
    scale_x_continuous(limits=c(0, l), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-top*.1, top), expand = c(0, 0)) +
    geom_vline(xintercept = c(1, l), size=.25) +
    geom_hline(yintercept = c(-top*.1, top), size=.25)

    grid.arrange(plot1, ncol=1)
  dev.off()

  #Plot histogram of bin counts
  jpeg(filename=paste(lab[k], "_counts", plot_sub, ".jpeg", sep=""), width=2500, height=1500)
    par(mar = c(7.0, 7.0, 7.0, 3.0))

    temp=sort(raw)[round(l*.01) : (l-round(l*.01))] 
    reads = hist(temp, breaks=100, plot=FALSE)
    # plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", xaxt="n", cex.main=3, cex.axis=2, cex.lab=2)
    plot(reads, col='black', main=paste("Frequency of Bin Counts for Sample ", lab[k], "\n(both tails trimmed 1%)", sep=""), xlab="Read Count (reads/bin)", cex.main=3, cex.axis=2, cex.lab=2)
    # axis(side=1, at=seq(min(temp), round(diff(range(temp))/20)*22, round(diff(range(temp))/20)), cex.axis=2)
    tu = par('usr')
    par(xpd=FALSE)
    clip(tu[1], as.integer(mean(temp)-(diff(reads$mids)/2)), tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(mean(temp)+(diff(reads$mids)/2), tu[2], tu[3], tu[4])
    plot(reads, col='gray50', add=TRUE)
    clip(tu[1], mean(temp) - sd(temp), tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(mean(temp) + sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray75', add=TRUE)
    clip(tu[1], mean(temp) - 2*sd(temp), tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE)
    clip(mean(temp) + 2*sd(temp), tu[2], tu[3], tu[4])
    plot(reads, col='gray90', add=TRUE) 
    legend("topright", inset=.05, legend=c("mean", "< 1σ", "> 1σ", "> 2σ"), fill=c("black", "gray50", "gray75", "gray90"), cex=2.5)
  dev.off()

  #Plot lorenz curves
  jpeg(filename=paste(lab[k], "_lorenz", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  nReads=sum(raw)
  uniq=unique(sort(raw))
  
  lorenz=matrix(0, nrow=length(uniq), ncol=2)
  a=c(length(which(raw==0)), tabulate(raw, nbins=max(raw)))
  b=a*(0:(length(a)-1))
  for (i in 2:length(uniq)) {
    lorenz[i,1]=sum(a[1:uniq[i]])/l
    lorenz[i,2]=sum(b[2:uniq[i]])/nReads
  }

  # smooth.spline needs >= 4 points...
  fit = data.frame(x=lorenz[,1], y=lorenz[,2])
  if(nrow(lorenz) >= 4)
  {
    spline = try(smooth.spline(lorenz))
    if(class(spline) != "try-error")
      fit = data.frame(x=spline$x, y=spline$y)
  }

  perf=data.frame(x=c(0,1), y=c(0,1))

  plot1 = try(ggplot() +
    geom_line(data=perf, aes(x=x, y=y, color="Perfect Uniformity"), size=3) +
    geom_line(data=fit, aes(x=x, y=y, color="Sample Uniformity"), size=3) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, .1)) +
    labs(title=paste("Lorenz Curve of Coverage Uniformity for Sample ", lab[k], sep=""), x="Cumulative Fraction of Genome", y="Cumulative Fraction of Total Reads") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    theme(legend.title=element_blank(), legend.text=element_text(size=40)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.15, .85)) +
    scale_color_manual(name='', values=c('Perfect Uniformity'="black", 'Sample Uniformity'=colors[cp,1])))

    grid.arrange(plot1, ncol=1)
  dev.off()

  ##Plot GC correction
  jpeg(filename=paste(lab[k], "_GC", plot_sub, ".jpeg", sep=""), width=2500, height=1250)

  low = lowess(GC[,1], log(normal2[,k]), f=0.05)
  app = approx(low$x, low$y, GC[,1])
  cor = exp(log(normal2[,k]) - app$y)
  
  uncorrected = data.frame(x=GC[,1], y=log(normal2[,k]))
  corrected = data.frame(x=GC[,1], y=log(cor))
  fit = data.frame(x=app$x, y=app$y)

  try(plot1 <- ggplot() +
    geom_point(data=uncorrected, aes(x=x, y=y), size=3) +
    geom_line(data=fit, aes(x=x, y=y, color="Lowess Fit"), size=3) +
    scale_x_continuous(limits=c(min(.3, min(GC[,1])), max(.6, max(GC[,1]))), breaks=seq(.3,.6,.05)) +
    labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Uncorrected)", sep=""), x="GC content", y="Normalized Read Counts (log scale)") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.85, .9)) +
    scale_color_manual(name='', values=colors[cp,1]))

  try(plot2 <- ggplot() +
    geom_point(data=corrected, aes(x=x, y=y), size=3) +
    scale_x_continuous(limits=c(min(.3, min(GC[,1])), max(.6, max(GC[,1]))), breaks=seq(.3,.6,.05)) +
    labs(title=paste("GC Content vs. Bin Count\nSample ", lab[k], " (Corrected)", sep=""), x="GC content", y="") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')))

    try(grid.arrange(plot1, plot2, ncol=2))
  dev.off()

  #Plot Scaled/Normalized Bin Count Histogram
  jpeg(filename=paste(lab[k], "_hist", plot_sub, ".jpeg", sep=""), width=2500, height=1500)

  clouds=data.frame(x=normal[,k]*CN)

  plot1 = ggplot() +
    geom_histogram(data=clouds, aes(x=x), binwidth=.05, color="black", fill="gray60") +
    geom_vline(xintercept=seq(0,10,1), size=1, linetype="dashed", color=colors[cp,1]) +
    scale_x_continuous(limits=c(0,10), breaks=seq(0,10,1)) +
    labs(title=paste("Frequency of Bin Counts for Sample \"", lab[k], "\"\nNormalized and Scaled by Predicted CN (", CNmult[1,k], ")", sep=""), x="Copy Number", y="Frequency") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black'))

    grid.arrange(plot1, ncol=1)
  dev.off()

  #Plot sum of squares error for each potential copy number
  jpeg(filename=paste(lab[k], "_SoS.jpeg", sep=""), width=2500, height=1500)

  top = max(outerColsums[,k])
  dat = data.frame(x=CNgrid, y=outerColsums[,k])
  lim = cbind(c(seq(0,5000,500), 1000000), c(50, 100, 100, 200, 250, 400, 500, 500, 600, 600, 750, 1000))
  step = lim[which(top<lim[,1])[1],]
  minSoS = data.frame(x=CNmult[1,k], y=CNerror[1,k])
  bestSoS = data.frame(x=CN, y=outerColsums[which(CNgrid==CN),k])

  plot1 = ggplot() +
    geom_line(data=dat, aes(x=x, y=y), size=3) +
    geom_point(data=dat, aes(x=x, y=y), shape=21, fill="black", size=5) +
    geom_point(data=minSoS, aes(x=x, y=y*1.02, color="Minimum SoS Error"), shape=18, size=15) +
    geom_point(data=bestSoS, aes(x=x, y=y*.98, color="Chosen Ploidy"), shape=18, size=15) +
    scale_x_continuous(limits=c(1.5, 6), breaks=seq(1.5, 6, .5)) +
    scale_y_continuous(limits=c(.5*min(outerColsums[,k]), top), breaks=seq(0, step[1], step[2])) +
    labs(title="Sum of Squares Error Across Potential Copy Number States", x="Copy Number Multiplier", y="Sum of Squares Error") +
    theme(plot.title=element_text(size=45, vjust=1.5)) +
    theme(axis.title.x=element_text(size=45, vjust=-2.8), axis.title.y=element_text(size=45, vjust=.1)) +
    theme(axis.text=element_text(color="black", size=45), axis.ticks=element_line(color="black")) +
    theme(plot.margin=unit(c(.5,1,1,1.5),"cm")) +
    theme(panel.background = element_rect(color = 'black')) +
    theme(legend.title=element_blank(), legend.text=element_text(size=45)) +
    theme(legend.key.height=unit(4,"line"), legend.key.width=unit(4,"line")) +
    theme(legend.position=c(.85, .9)) +
    scale_color_manual(name='', values=c('Minimum SoS Error'=colors[cp,1], 'Chosen Ploidy'=colors[cp,2]))

    grid.arrange(plot1, ncol=1)
  dev.off()

  t5 <- Sys.time()
  total_plots <<- total_plots + (t5-t4)

  v <- list("out" = outliers, "breaks" = breaks, "loc" = loc, "raw" = raw, "normal" = normal, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
  return(v)
  #return(outliers)
}

write_results <- function(all_loc, all_bounds, all_normal, all_fixed, all_final, all_breaks, stats, desc='') {
  sink("results.txt", append=TRUE)
  cat(paste("Writing results", "\n", sep=''))
  # Store processed sample information
  loc2=all_loc
  loc2[,3]=loc2[,2]
  pos = cbind(c(1,all_bounds[,2]), c(all_bounds[,2], l))

  #
  for (i in 1:nrow(pos))
  {
    cat(paste(paste(pos[i,], collapse=","), "\n", sep=''))
    # If only 1 bin in a chromosome
    if( (pos[i,2] - pos[i,1]) == 0 ) {
      loc2[pos[i,1],1] = 1
    # If two bins.......
    } else if( (pos[i,2] - pos[i,1]) == 1 ) {
      loc2[pos[i,1],1] = 1
      loc2[pos[i,2],1] = loc2[pos[i,1],2] + 1
    } else {
      loc2[pos[i,1]:(pos[i,2]-1),2]=c(1,loc[pos[i,1]:(pos[i,2]-2),2]+1)
    }
  }

  # 
  loc2[nrow(loc2),2]=loc2[nrow(loc2)-1,3]+1
  colnames(loc2)=c("CHR","START", "END")

  suffix <- ""
  if (clusterCells) {
    suffix <- "Clust"
  }
  
  cat(paste("Writing", "\n", sep=''))
  write.table(cbind(loc2,all_normal), file=paste(user_dir, "/SegNorm", suffix,  sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
  write.table(cbind(loc2,all_fixed), file=paste(user_dir, "/SegFixed", suffix, sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
  write.table(cbind(loc2,all_final), file=paste(user_dir, "/SegCopy", suffix, sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
  write.table(cbind(loc2,all_breaks), file=paste(user_dir, "/SegBreaks", suffix, sep=""), row.names=FALSE, col.names=c(colnames(loc2),lab), sep="\t", quote=FALSE)
  write.table(stats, file=paste(user_dir, "/SegStats", suffix, sep=""), sep="\t", quote=FALSE)
  cat(paste("Done", "\n", sep=''))
  sink()
}

compare_CN <- function(all_final, cell1, cell2, cell_ids, thresh=30) {
    num_mismatches = 0
    CN1 <- all_final[,cell_ids[[cell1]]]
    CN2 <- all_final[,cell_ids[[cell2]]]
    num_bins = length(CN1)
    curr_mismatch_len = 0
    for (i in 1:num_bins) {
        if (CN1[i] == CN2[i]) {
            # Reset segment counter
            curr_mismatch_len = 0
        } else {
            curr_mismatch_len = curr_mismatch_len + 1
            if (curr_mismatch_len == thresh)
                num_mismatches = num_mismatches + 1
        }
    }
    return(num_mismatches)
}

score_cluster <- function(all_final, clust, threshold, cell_ids) {
  # Count the number of large mismatches among all pairs of same-cluster cells when clustered with the given threshold
  memb <- cutree(clust, h=threshold)
  groups <- tapply(names(memb), memb, c)

  #cat(threshold, "\n")

  mismatches <- 0
  for (n in names(groups)) {
    curr_clust = groups[[n]]
    clust_size = length(curr_clust)
    if (clust_size < 2)
      next

    #cat(paste("  ", clust_size, "\n",sep=''))
    for (i in 1:(clust_size-1)) {
      for (j in (i+1):clust_size) {
        mismatches = mismatches + compare_CN(all_final, curr_clust[i], curr_clust[j], cell_ids)
      }
    }
  }
  return(mismatches)
}

gen_clusts <- function(all_fixed, all_final, desc='') {
  # Calculate read distance matrix for clustering
  d = dist(t(all_fixed), method=dm)
  if(cm == "NJ")
  {
    library(ape)
    clust = nj(d)
    clust$tip.label = lab
    if (!clusterCells) {
      saveRDS(clust, file=paste(user_dir, "/clust.hclust", sep=""))
    }
    write.tree(clust, file=paste(user_dir, "/clust", desc, outSuffix, ".newick", sep=""))
  } else {
    clust = hclust(d, method = cm)
    clust$labels = lab  
    if (!clusterCells) {
      saveRDS(clust, file=paste(user_dir, "/clust.hclust", sep=""))
    }
    write(hc2NewickFixed(clust), file=paste(user_dir, "/clust", desc, outSuffix, ".newick", sep=""))
  }

  cell_ids <- list()
  for(k in 1:w)
  {
    cell_ids[[lab[k]]] <- k
  }

  if (!clusterCells) {
    #for (i in seq(3,20,by=.1)) {
    #  memb <- cutree(clust, h=i)
    #  groups <- tapply(names(memb), memb, c)
    #  cat(i, ":\t", length(names(groups)), "\n")
    #}

    # Predict optimal clustering threshold
    maxThresh = 0
    maxScore = 0
    # Find thresholds s.t. maxThresh = minThresh+1, score(minThresh) = 0 and score(maxThresh) > 0
    while (maxScore == 0) {
      minThresh = maxThresh
      minScore = 0
      maxThresh = maxThresh + 1
      maxScore = score_cluster(all_final, clust, maxThresh, cell_ids)
      if (maxThresh > 20)
        break
    }

    step = .1
    for (i in seq(minThresh+step,maxThresh-step,by=step)) {
      score = score_cluster(all_final, clust, i, cell_ids)
      if (score == 0) {
        minThresh = i
        minScore = 0
      } else {
        maxThresh = i
        maxScore = score
        break
      }
    }

    write(minThresh, file="optimal_thresh.txt")
  }

  ###
  command=paste("java -cp ", main_dir, "/forester_1025.jar org.forester.application.phyloxml_converter -f=nn ", user_dir, "/clust", desc, outSuffix, ".newick ", user_dir, "/clust", outSuffix, ".xml", sep="");
  unlink( paste(user_dir, "/clust", outSuffix, ".xml", sep="") );
  system(command);
  ###

  #Plot read cluster
  jpeg(paste("clust", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
    op = par(bg = "gray85")
    plot(clust, xlab="Sample", hang=-1, ylab=paste("Distance (", dm, ")", sep=""), lwd=2, cex.main=2, cex.lab=2, cex.axis=2)
  dev.off()

  pdf(paste("clust", outSuffix, ".pdf", sep=""), width=10, height=7)
    op = par(bg = "gray85")
    plot(clust, xlab="Sample", ylab=paste("Distance (", dm, ")", sep=""), lwd=2)
  dev.off()

  pdf(paste("clust_aligned", outSuffix, ".pdf", sep=""), width=10, height=7)
    op = par(bg = "gray85")
    plot(clust, xlab="Sample", hang=-1, ylab=paste("Distance (", dm, ")", sep=""), lwd=2)
  dev.off()

  statusFile=file( paste(user_dir, "/", status, sep="") )
  writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>Computing Cluster (Copy Number)</processingfile>", sep=""), paste("<percentdone>", ((w+2)*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
  close(statusFile)

  #Calculate copy number distance matrix for clustering
  d2 = dist(t(all_final), method = dm)
  #clust2 = hclust(d2, method = cm)
  #clust2$labels = lab
  #write(hc2NewickFixed(clust2), file=paste(user_dir, "/clust2.newick", sep=""))
  if(cm == "NJ"){
    library(ape)
    clust2 = nj(d2)
    clust2$tip.label = lab
    write.tree(clust2, file=paste(user_dir, "/clust2", outSuffix, ".newick", sep=""))
  }else{
    clust2 = hclust(d2, method = cm)
    clust2$labels = lab  
    write(hc2NewickFixed(clust2), file=paste(user_dir, "/clust2", outSuffix, ".newick", sep=""))
  }

  ###
  command=paste("java -cp ", main_dir, "/forester_1025.jar org.forester.application.phyloxml_converter -f=nn ", user_dir, "/clust2", outSuffix, ".newick ", user_dir, "/clust2", outSuffix, ".xml", sep="");
  unlink( paste(user_dir, "/clust2", outSuffix, ".xml", sep="") );
  system(command);
  ###

  #Plot copy number cluster
  jpeg(paste("clust2", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
    op = par(bg = "gray85")
    plot(clust2, xlab="Sample", hang=-1, ylab=paste("Distance (", dm, ")", sep=""), lwd=2)
  dev.off()

  pdf(paste("clust2", outSuffix, ".pdf", sep=""), width=10, height=7)
    op = par(bg = "gray85")
    plot(clust2, xlab="Sample", hang=-1, ylab=paste("Distance (", dm, ")", sep=""), lwd=2)
  dev.off()

  #Calculate correlation distance matrix for clustering
  d3 = as.dist((1 - cor(all_final))/2)
  #clust3=hclust(d3, method = cm)
  #clust3$labels=lab
  #write(hc2NewickFixed(clust3), file=paste(user_dir, "/clust3.newick", sep=""))
  if(cm == "NJ"){
    library(ape)
    clust3 = nj(d3)
    clust3$tip.label = lab
    write.tree(clust3, file=paste(user_dir, "/clust3", outSuffix, ".newick", sep=""))
  }else{
    clust3 = hclust(d3, method = cm)
    clust3$labels = lab  
    write(hc2NewickFixed(clust3), file=paste(user_dir, "/clust3", outSuffix, ".newick", sep=""))
  }

  ###
  command=paste("java -cp ", main_dir, "/forester_1025.jar org.forester.application.phyloxml_converter -f=nn ", user_dir, "/clust3", outSuffix, ".newick ", user_dir, "/clust3", outSuffix, ".xml", sep="");
  unlink( paste(user_dir, "/clust3", outSuffix, ".xml", sep="") );
  system(command);
  ### 

  #Plot correlation cluster
  jpeg(paste("clust3", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
    op = par(bg = "gray85")
    plot(clust3, xlab="Sample", hang=-1, ylab="Distance (Pearson correlation)", lwd=2)
  dev.off()

  pdf(paste("clust3", outSuffix, ".pdf", sep=""), width=10, height=7)
    op = par(bg = "gray85")
    plot(clust3, xlab="Sample", hang=-1, ylab="Distance (Pearson correlation)", lwd=2)
  dev.off()

  clusts <- list("clust" = clust, "clust2" = clust2, "clust3" = clust3)
  return(clusts)
}

plot_heatmaps <- function(rawBPs, fixedBPs, finalBPs, clust, clust2, clust3) {
  colnames(rawBPs) = lab
  colnames(fixedBPs) = lab
  colnames(finalBPs) = lab

  # Need to root NJ tree, make tree ultrametric by extending branch lengths then convert to hclust object!
  phylo2hclust = function(phy)
  {
    # Root tree
    clustR = root(phy, outgroup=1, resolve.root=TRUE)
    # If edge lengths are exactly 0, chronopl will delete those edges.....
    clustRE= clustR$edge.length
    clustRE[which(clustRE == 0)] = clustRE[which(clustRE == 0)]+ 1e-4
    clustRE[which(clustRE < 0)] = 1e-4
    clustR$edge.length = clustRE
    # Chronopl to make tree ultrametric
    clustU = chronopl(clustR, 0)
    phy  = as.hclust(clustU)
    return(phy)
  }

  #
  if(cm == "NJ"){
    clust  = phylo2hclust(clust)
    clust2 = phylo2hclust(clust2)
    clust3 = phylo2hclust(clust3)
  }

  write("Making heatRaw.jpeg", stderr())
  jpeg(paste("heatRaw", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
  heatmap.2(t(rawBPs), Colv=FALSE, Rowv=as.dendrogram(clust), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=bluered(2))
  dev.off()

  write("Making heatNorm.jpeg", stderr())
  step=quantile(fixedBPs, c(.98))[[1]]
  jpeg(paste("heatNorm", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
  heatmap.2(t(fixedBPs), Colv=FALSE, Rowv=as.dendrogram(clust), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=bluered(15), breaks=seq(0,step,step/15))
  dev.off()

  write("Making heatCN.jpeg", stderr())
  step=min(20, quantile(finalBPs, c(.98))[[1]])
  jpeg(paste("heatCN", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
  heatmap.2(t(finalBPs), Colv=FALSE, Rowv=as.dendrogram(clust2), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=colorRampPalette(c("white","green","green4","violet","purple"))(15), breaks=seq(0,step,step/15))
  dev.off()

  write("Making heatCor.jpeg", stderr())
  jpeg(paste("heatCor", outSuffix, ".jpeg", sep=""), width=2000, height=1400)
  heatmap.2(t(finalBPs), Colv=FALSE, Rowv=as.dendrogram(clust3), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=colorRampPalette(c("white","steelblue1","steelblue4","orange","sienna3"))(15), breaks=seq(0,step,step/15))
  dev.off()

  #jpeg("heatRaw.jpeg", width=2000, height=1400)
  #heatmap3(t(rawBPs), distfun = function(x) dist(x), Colv=FALSE, Rowv=as.dendrogram(clust), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=bluered(2))
  #dev.off()
  #
  #step=quantile(fixedBPs, c(.98))[[1]]
  #jpeg("heatNorm.jpeg", width=2000, height=1400)
  #heatmap3(t(fixedBPs), distfun = function(x) dist(x), Colv=FALSE, Rowv=as.dendrogram(clust), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=bluered(15), breaks=seq(0,step,step/15))
  #dev.off()
  #
  #step=min(20, quantile(finalBPs, c(.98))[[1]])
  #jpeg("heatCN.jpeg", width=2000, height=1400)
  #heatmap3(t(finalBPs), distfun = function(x) dist(x), Colv=FALSE, Rowv=as.dendrogram(clust2), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=colorRampPalette(c("white","green","green4","violet","purple"))(15), breaks=seq(0,step,step/15))
  #dev.off()
  #
  #jpeg("heatCor.jpeg", width=2000, height=1400)
  #heatmap3(t(finalBPs), distfun = function(x) dist(x), Colv=FALSE, Rowv=as.dendrogram(clust3), margins=c(5,20), dendrogram="row", trace="none", xlab="Bins", ylab="Samples", cex.main=2, cex.axis=1.5, cex.lab=1.5, cexCol=.001, col=colorRampPalette(c("white","steelblue1","steelblue4","orange","sienna3"))(15), breaks=seq(0,step,step/15))
  #dev.off()
}


#splitBins <- FALSE
num_iters = 4
all_loc <- loc
all_raw <- raw
all_normal <- normal

all_fixed <- matrix(0,l,w)
all_final <- matrix(0,l,w)
all_breaks <- normal
all_bounds <- bounds

for(k in 1:w)
{
  cat('===',k,'===\n')
  if (!improveBounds) {
    v = processSample(k, loc, raw[,k], GC[,1], bounds, plot_sub='')
  } else {
    if (split_mode == 3) {
      v = processSampleRefine(k, loc, curr_raw, GC[,1], bounds, plot_sub='')
    } else {
      v = processSample(k, loc, curr_raw, GC[,1], bounds, plot_sub='')
      for (iter in 1:num_iters) {
        cat(paste("Iter ", iter, "\n", sep=''))
        outliers <- v$out
        breaks <- v$breaks
        #v = processSample(k, v$loc, v$raw, v$gc, v$bounds, plot_sub=paste(".v", iter+1, sep=''), outliers, v$segbreaks, v$wgt)
        v = processSampleNoNewBins(k, loc, curr_raw, GC[,1], bounds, outliers, breaks, iter=iter+1, plot_sub=paste(".v", iter+1, sep=''))
      }
    }
  }

  # TODO: Add bounds(?)
  curr_loc <- v$loc
  curr_raw <- v$raw
  curr_normal <- v$normal
  curr_fixed <- v$fixed
  curr_final <- v$final
  curr_breaks <- v$segbreaks

  #cat("Curr normal:\n")
  #cat(paste(paste(curr_normal, sep=''), "\n", sep=''))
  #cat("Curr fixed:\n")
  #cat(paste(length(curr_fixed), " x ", length(curr_fixed[,1]), "\n", sep=''))
  #cat("All fixed:\n")
  #cat(paste(length(all_fixed), " x ", length(all_fixed[,1]), "\n", sep=''))

  #num_bins <- length(all_loc)
  #i <- 1
  #while (i <= num_bins) {
  #  if (all_loc[i,1] != curr_loc[i,1]) {
  #    cat("Error! Tried to split bin at chromsome boundary!\n")
  #    quit()
  #  }
  #  if (all_loc[i,2] > curr_loc[i,2]) {
  #    # Split existing bin
  #    all_loc <- append(all_loc, curr_loc[i,], after=i-1)
  #    # Read counts are weighted by bin width, so splitting a bin in half does not change the value of either bin
  #    all_raw <- append(all_raw, all_raw[i,], after=i-1)
  #    all_normal <- append(all_normal, all_normal[i,], after=i-1)
  #    all_fixed <- append(all_fixed, all_fixed[i,], after=i-1)
  #    all_final <- append(all_final, all_final[i,], after=i-1)
  #    all_breaks <- append(all_breaks, all_breaks[i,], after=i-1)
  #    all_breaks[i+1,] <- matrix(0,1,w)

  #    for (chrom_id in 1:length(all_bounds[,1]))
  #      if (all_bounds[chrom_id,2] >= i)
  #        all_bounds[chrom_id,2] = all_bounds[chrom_id,2]+1

  #    num_bins <- num_bins+1      
  #  }
  #  else if (all_loc[i,2] < curr_loc[i,2]) {
  #    # Split new bin
  #    curr_loc <- append(curr_loc, all_loc[i,], after=i-1)
  #    curr_raw <- append(curr_raw, curr_raw[i], after=i-1)
  #    curr_normal <- append(curr_normal, curr_normal[i], after=i-1)
  #    curr_fixed <- append(curr_fixed, curr_fixed[i], after=i-1)
  #    curr_final <- append(curr_final, curr_final[i], after=i-1)
  #    curr_breaks <- append(curr_breaks, curr_breaks[i], after=i-1)
  #    curr_breaks[i+1] <- 0
  #  }
  #  i <- i+1
  #}

  all_raw[,k] <- curr_raw
  all_normal[,k] <- curr_normal
  #all_normal[,k] <- curr_normal[,k]
  all_fixed[,k] <- curr_fixed
  all_final[,k] <- curr_final
  all_breaks[,k] <- curr_breaks
}

sink("timing.txt", append=TRUE)
timeD <- proc.time()
cat(paste((timeD-timeC)['elapsed'], "\n", sep=""))
cat(paste("  ", total_seg, "\n", sep=""))
cat(paste("  ", total_plots, "\n", sep=""))
cat(paste("  ", total_other, "\n", sep=""))
sink("results.txt", append=TRUE)

# ------------------------------------------------------------------------------
# -- Save processed data
# ------------------------------------------------------------------------------

# Update status
statusFile=file( paste(user_dir, "/", status, sep="") )
writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>Saving Data</processingfile>", sep=""), paste("<percentdone>", (w*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
close(statusFile)

# Close output stream
sink()

write_results(all_loc, all_bounds, all_normal, all_fixed, all_final, all_breaks, stats)
sink("timing.txt", append=TRUE)
timeE <- proc.time()
cat(paste((timeE-timeD)['elapsed'], "\n", sep=""))
sink("results.txt", append=TRUE)

# ------------------------------------------------------------------------------
# -- Generate phylo trees
# ------------------------------------------------------------------------------

statusFile = file( paste(user_dir, "/", status, sep="") )
writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>Computing Cluster (Read Count)</processingfile>", sep=""), paste("<percentdone>", ((w+1)*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
close(statusFile)

# Ignore sex chromosomes if specified
# TODO: Recompute l for added split bins
if (sex == 0) {
  l=all_bounds[(dim(all_bounds)-1)[1],][[2]]-1
  all_raw=all_raw[1:l,]
  all_final=all_final[1:l,]
  all_fixed=all_fixed[1:l,]
  all_breaks=all_breaks[1:l,]
  all_normal=all_normal[1:l,]
}

hc2NewickFixed <- function(hc, flat=TRUE) {
  distA <- 0
  distB <- 0
  if (is.null(hc$labels))
    labels <- seq(along=hc$order)
  else
    labels <- hc$labels

  putparenthesis <- function(i) {
    ## recursive function
    ## i: index of row in hc$merge
    j <- hc$merge[i, 1]
    k <- hc$merge[i, 2]
    
    if (j < 0) {
      left <- labels[-j]
      distA <- hc$height[i]
    } else {
      left <- putparenthesis(j)
      distA <- hc$height[i] - hc$height[j]
    }

    if (k < 0) {
      right <- labels[-k]
      distB <- hc$height[i]
    } else {
      right <- putparenthesis(k)
      distB <- hc$height[i] - hc$height[k]
    }

    if (flat)
      return(paste("(", left, ":", distA, ",", right, ":", distB, ")", sep=""))
    else
      return(list(left=left, right=right, dist=dist))
  }
  
  n <- putparenthesis(nrow(hc$merge))
  if (flat)
    n <- paste(n, ";", sep="")
  
  return(n)
}

cat(paste("Generating clusters", "\n", sep=''))
clusts <- gen_clusts(all_fixed, all_final, '')
sink("timing.txt", , append=TRUE)
timeF <- proc.time()
cat(paste((timeF-timeE)['elapsed'], "\n", sep=""))
sink("results.txt", append=TRUE)

# ------------------------------------------------------------------------------
# -- Generate heatmaps
# ------------------------------------------------------------------------------

statusFile=file( paste(user_dir, "/", status, sep="") )
writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>Creating Heat Maps</processingfile>", sep=""), paste("<percentdone>", ((w+3)*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
close(statusFile)

plot_heatmaps(all_breaks, all_fixed, all_final, clusts$clust, clusts$clust2, clusts$clust3)
sink("timing.txt", append=TRUE)
timeG <- proc.time()
cat(paste((timeG-timeF)['elapsed'], "\n", sep=""))
sink("results.txt", append=TRUE)

statusFile=file( paste(user_dir, "/", status, sep="") )
writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>Finished</processingfile>", sep=""), paste("<percentdone>100</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
close(statusFile)

#sink("timing.txt", append=TRUE)
#timeH <- Sys.time()
#cat(paste(timeH-timeA, "\n", sep=""))
