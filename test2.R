processSample <- function(k, loc, raw, GC, bounds, plot_sub=NULL, outliers=NULL, breaks=NULL, bin_wgts=NULL) {
  statusFile = file( paste(user_dir, "/", status, sep="") )
  writeLines(c("<?xml version='1.0'?>", "<status>", "<step>3</step>", paste("<processingfile>", lab[k], "</processingfile>", sep=""), paste("<percentdone>", (k*100)%/%(w+4), "</percentdone>", sep=""), "<tree>clust.xml</tree>", "</status>"), statusFile)
  close(statusFile)

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
        bin_id = outliers[i]
        if (loc[bin_id+i-2,1] != loc[bin_id+i-1,1])
          next 

        bin_start = loc[bin_id+i-2,2]
        bin_end = loc[bin_id+i-1,2]
        bin_mid = (bin_end - bin_start) * breaks[i] + bin_start

        fine_start = 1
        #while (loc_fine[fine_start,1] != loc[bin_id+i-1,1]) # Find the right chromosome
        #  fine_start = fine_start+1
        while (loc_fine[fine_start,1] != loc[bin_id+i-1,1] || loc_fine[fine_start+1,2] <= bin_start) # Find the right bin
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
        wgt_left = bins_left / (bins_left+bins_right)
        wgt_right = bins_right / (bins_left+bins_right)

        curr_wgt = bin_wgts[bin_id+i-1]
        bin_wgts[bin_id+i-1] = curr_wgt/wgt_left
        bin_wgts = append(bin_wgts, curr_wgt/wgt_right, after=bin_id+i-1)
        #loc = append(loc, loc_fine[fine_mid,], after=bin_id+i-2)
        loc = rbind(loc[1:(bin_id+i-2),], loc_fine[fine_mid,], loc[(bin_id+i-1):length(loc[,1]),])

        bin_raw = raw[bin_id+i-1]

        left_raw = sum(raw_fine[(fine_start+2):fine_mid,k]) + raw_fine[(fine_start+1),k] * left_bin_overlap
        right_raw = sum(raw_fine[(fine_mid+1):(fine_end-1),k]) + raw_fine[fine_end,k] * right_bin_overlap
        raw[bin_id+i-1] = left_raw

        raw = append(raw, right_raw, after=bin_id+i-1)

        #GC = rbind(GC[1:(bin_id+i-1),], GC[(bin_id+i-1),], GC[(bin_id+i-1):length(GC[,1]),])
        GC = append(GC, GC[(bin_id+i-1)], after=bin_id+i-1)
        F = append(F, F[(bin_id+i-1)], after=bin_id+i-1)
        for (chrom_id in 1:length(bounds[,1]))
          if (bounds[chrom_id,2] > bin_id)
            bounds[chrom_id,2] = bounds[chrom_id,2]+1
      }
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

  # Determine breakpoints and extract chrom/locations
  CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
  CNA.smoothed = smooth.CNA(CNA.object)
  segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
  frag         = segs$output[,2:3]

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
    fixed[bps[i]:(bps[i+1]-1)] = curr_cn

    cn_change = curr_cn - last_cn

    num_b = length(bounds[,1])
    while (chrom_i <= num_b && bounds[chrom_i,2] < bps[i])
      chrom_i = chrom_i + 1

    if ((chrom_i <= num_b) && (bounds[chrom_i,2] != bps[i]) && (abs(cn_change) >= cn_change_threshold)) {
      z1 = (normal[bps[i]-1] - last_mean) / last_sd
      z2 = (normal[bps[i]] - curr_mean) / curr_sd
      found_outlier = F
      if (curr_cn > last_cn) {
        if ((z1 > zthresh) && (z2 > 0 || z1 > -z2)) {
          outliers[j] = bps[i]-1
          found_outlier = T
        }
        else if (z2 < -zthresh) {
          outliers[j] = bps[i]
          found_outlier = T
        }
      }
      else {
        if ((z1 < -zthresh) && (z2 < 0 || -z1 > z2)) {
          outliers[j] = bps[i]-1
          found_outlier = T
        }
        else if (z2 > zthresh) {
          outliers[j] = bps[i]
          found_outlier = T
        }
      }
      if (found_outlier) {
        split_mode = 1
        if (split_mode == 0) {
          # Split in half
          breaks[j] = 0.5
        }
        else if (split_mode == 1) {
          # Simple interpolation
          #breaks[j] = (normal[outliers[j]] - last_mean) / (curr_mean - last_mean)
          breaks[j] = (curr_mean - normal[outliers[j]]) / (curr_mean - last_mean)
        }
        else if (split_mode == 2) {
          # z-score-based interpolation
          zscore1 = abs((normal[outliers[j]] - last_mean)) / last_sd
          zscore2 = abs((normal[outliers[j]] - curr_mean)) / curr_sd
          breaks[j] = zscore2 / (zscore1 + zscore2)
        }
        j = j+1
      }
    }

    last_cn = curr_cn
    last_mean = curr_mean
    last_sd = curr_sd
  }
  #cat(paste("\nOutliers: ", paste(outliers, collapse=","), "\n", sep=''))
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


  write.table(cbind(loc,final), file=paste(lab[k], "_CN", plot_sub, ".tsv", sep=""), row.names=FALSE, col.names=c(colnames(loc),lab[k]), sep="\t", quote=FALSE)

  #if (T) {
  if (missing(plot_sub)) {
    if (length(raw) == length(bin_wgts)) {
      for (raw_i in 1:length(raw))
        raw[raw_i] = raw[raw_i] / bin_wgts[raw_i]
    }

    v <- list("out" = outliers, "breaks" = breaks, "wgt" = bin_wgts, "loc" = loc, "gc" = GC, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
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
  clouds=data.frame(x=1:length(normal), y=normal*CN)
  amp=data.frame(x=which(final>2), y=final[which(final>2)])
  del=data.frame(x=which(final<2), y=final[which(final<2)])
  flat=data.frame(x=which(final==2), y=final[which(final==2)])
  anno=data.frame(x=(pos[,2]+pos[,1])/2, y=-top*.05, chrom=substring(c(as.character(bounds[,1]), "chrY"), 4 ,5))
  # TODO: Fix outliers, array contains NAs and is sometimes NULL
  print("Outliers:")
  print(outliers)
  outlier_pts=data.frame(x=outliers, y=normal[outliers]*CN)

  plot_chroms <- FALSE
  if (plot_chroms) {
    for (chr_id in 1:(length(pos[,1])-1)) {
      #cat(paste(pos[chr_id], ",", pos[chr_id+1], "\n", sep=''))
      jpeg(filename=paste(lab[k], "_CN", plot_sub, ".chr", chr_id, ".jpeg", sep=""), width=3000, height=750)
      plot1 = ggplot() +
        geom_rect(data=rectangles1, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray85', alpha=0.75) +
        geom_rect(data=rectangles2, aes(xmin=X1, xmax=X2, ymin=-top*.1, ymax=top), fill='gray75', alpha=0.75) +
        geom_point(data=clouds, aes(x=x, y=y), color='gray45', size=3) +
        geom_point(data=flat, aes(x=x, y=y), size=4) +
        geom_point(data=amp, aes(x=x, y=y), size=4, color=colors[cp,1]) +
        geom_point(data=del, aes(x=x, y=y), size=4, color=colors[cp,2]) +
        geom_point(data=outlier_pts, aes(x=x, y=y), size=4, color='blue') +
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
        geom_hline(yintercept = c(-top*.1, top), size=.5)+
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
      #geom_point(data=outlier_pts, aes(x=x, y=y), size=4, color='blue') +
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
  #jpeg(filename=paste(lab[k], "_GC", plot_sub, ".jpeg", sep=""), width=2500, height=1250)

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

  if (length(raw) == length(bin_wgts)) {
    for (raw_i in 1:length(raw))
      raw[raw_i] = raw[raw_i] / bin_wgts[raw_i]
  }
  v <- list("out" = outliers, "breaks" = breaks, "wgt" = bin_wgts, "loc" = loc, "gc" = GC, "raw" = raw, "normal" = normal*CN, "fixed" = fixed, "final" = final, "bounds" = bounds, "segbreaks"=seg_breaks)
  return(v)
  #return(outliers)
}
