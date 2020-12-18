# Script for tree swallow RRBS data from Ithaca NY
# Written by Conor Taff: cct663@gmail.com ~ cct63@cornell.edu
# Last updated 11/2020

# Largely based on vignette here: http://127.0.0.1:18146/library/methylKit/doc/methylKit.html

# Load libraries ----
    pacman::p_load("methylKit", "tidyverse", "here", "gridExtra", "ggpubr", "viridis", "DSS", "genomation",
                   "lme4", "lmerTest", "reshape2")

# Load data ----
    cov_list <- list.files(here("0_processed_data/bismark_cov_output"))
    d_bis <- read.delim(here("0_processed_data/bismark_sample_summary.txt"))
    d_sample <- read.delim(here("1_raw_data/rrbs_sample_metadata.txt"))
    
# Set colors ----
    cort_col <- "#56B4E9"
    cont_col <- "#E69F00"
    pre_col <- "#999999"
    
# Simple plots of reads and methylation ----
    # reads and methylation levels
            p_a <- d_bis %>%
                pivot_longer(cols = c("Total.Reads", "Aligned.Reads"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, binwidth = 1000, position = "identity", color = "gray30") +
                    theme_classic() + xlab("Total Sequences / 1000") +
                    scale_fill_manual(values = c(Total.Reads = "#F0E442", Aligned.Reads = "#56B4E9"),
                                      labels = c("Aligned Reads", "Total Reads")) +
                    ylab("Number of Samples") +
                    theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "A")
            
            p_d <- d_bis %>%
                pivot_longer(cols = c("Methylated.CpGs", "Unmethylated.CpGs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, binwidth = 1000, position = "identity", color = "gray30") +
                    theme_classic() + xlab("Number of CpGs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CpGs = "#009E73", Methylated.CpGs = "#E69F00"),
                                      labels = c("Methylated", "Unmethylated")) +
                    ylab("Number of Samples") +
                    theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "D")
            
            p_b <- d_bis %>%
                pivot_longer(cols = c("Methylated.CpHs", "Unmethylated.CpHs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, binwidth = 1000, position = "identity", color = "gray30") +
                    theme_classic() + xlab("Number of CpHs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CpHs = "#009E73", Methylated.CpHs = "#E69F00"),
                                      labels = c("Methylated", "Unmethylated")) +
                    ylab("Number of Samples") + ylim(c(0, 130)) +
                    theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "B") +
                    annotate("text", x = 38e3, y = 25, label = "1.0% methylation", size = 3)
            
            p_c <- d_bis %>%
                pivot_longer(cols = c("Methylated.CHHs", "Unmethylated.CHHs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, binwidth = 1000, position = "identity", color = "gray30") +
                    theme_classic() + xlab("Number of CHHs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CHHs = "#009E73", Methylated.CHHs = "#E69F00"),
                                      labels = c("Methylated", "Unmethylated")) +
                    ylab("Number of Samples") + ylim(c(0, 120)) +
                    theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "C") +
                    annotate("text", x = 75e3, y = 25, label = "0.6% methylation", size = 3)
            
            p_e <- ggplot(data = d_bis, mapping = aes(x = (Methylated.CpGs / (Unmethylated.CpGs + Methylated.CpGs)) * 100)) +
                geom_histogram(fill = "#D55E00", binwidth = 2.5, color = "gray30", alpha = 0.6) + 
                theme_classic() + xlab("CpG Methylation Percentage") +
                ylab("Number of Samples") +
                annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "E")
            
            ggsave(here("3_markdown_summary/read_summary.png"),
                grid.arrange(p_a, p_b, p_c, p_d, p_e,
                         layout_matrix = rbind(c(1, 1, 2, 4),
                                               c(1, 1, 3, 5))),
                device = "png", width = 12.1, height = 5.9)
    
    
    # Reads vs. alignment
        p_a <- ggplot(data = d_bis, mapping = aes(x = Total.Reads / 1000, y = Aligned.Reads / 1000)) + 
            geom_point(col = "slateblue") + theme_classic() + 
            geom_smooth(method = "lm", col = "coral3") +
            xlab("Total Reads / 1000") + ylab("Aligned Reads / 1000") +
            annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "A") +
            annotate("text", x = 15e3, y = 1200, label = "51.1% Alignment", size = 3)
        
        p_b <- ggplot(data = d_bis, mapping = aes(x = Aligned.Reads / 1000, y = (Methylated.CpGs / (Unmethylated.CpGs + Methylated.CpGs)) * 100)) +
            geom_point(col = "slateblue") + theme_classic() + geom_smooth(method = "lm", col = "coral3") +
            annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5 ,label = "B") +
            xlab("Aligned Reads / 1000") + ylab("CpG Percent Methylation")
        
        ggsave(here("3_markdown_summary/reads_align_meth.png"), 
               ggarrange(p_a, p_b),
               device = "png", width = 9, height = 4.5)
            
      
    # Methylation conversion from spikes
        p1 <- ggplot(data = d_sample, mapping = aes(x = meth_conv)) +
            geom_histogram(fill = "gray70", binwidth = 0.25, col = "gray30") +
            theme_classic() + xlab("Methylation Conversion Percent") +
            ylab("Number of Samples") + ggtitle("Unmethylated Spike In") +
            geom_vline(xintercept = 2, linetype = "dotted", col = "coral3", size = 1.5) +
            geom_vline(xintercept = mean(d_sample$meth_conv), col = "slateblue", size = 1.5) +
            annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "A")
        
        p2 <- ggplot(data = d_sample, mapping = aes(x = unmeth_conv)) +
            geom_histogram(fill = "gray70", binwidth = 0.25, col = "gray30") +
            theme_classic() + xlab("Methylation Conversion Percent") +
            ylab("Number of Samples") + ggtitle("Methylated Spike In") +
            geom_vline(xintercept = 98, linetype = "dotted", col = "coral3", size = 1.5) +
            geom_vline(xintercept = mean(d_sample$unmeth_conv), col = "slateblue", size = 1.5) +
            annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "B")
        
        ggsave(here("3_markdown_summary/meth_conversion.png"),
               ggarrange(p1, p2),
               device = "png", width = 8, height = 4)
    
# Pull sample id ----
    # merge to file that has filename and sample id
        cov_l <- as.data.frame(cov_list)
        colnames(cov_l) <- "bis_processed_file"
        cov_l <- plyr::join(cov_l, d_sample, "bis_processed_file")
        
    # Create group labels
        d_sample$group_id <- paste(d_sample$group, d_sample$treatment, sep = "_")
        grps <- data.frame(group_id = c("1pre_cort", "1pre_control", "2post_cort", "2post_control", "3cross_control", "3cross_cort", "_"),
                           comp_grp = c("D", "A", "E", "B", "C", "F", "G"))
        d_sample <- plyr::join(d_sample, grps, "group_id")
        labels_t <- data.frame(treatment = c("cort", "control"), trt_label = c("Corticosterone", "Control"))
        d_sample <- plyr::join(d_sample, labels_t, "treatment")
        
        labels2 <- data.frame(comp_grp = c("A", "B", "C", "D", "E", "F"), cross_grp = c("CAD", "B", "CAD", "CAD", "E", "F"))
        d_sample <- plyr::join(d_sample, labels2, "comp_grp")

# Figure describing samples ----
     
    stage_names <- c(
      '1pre' = "Pre-Treatment",
      '2post' = "Post-Treatment",
      '3cross' = "Following Year")
        
    treat_names <- c(
      'control' = "Control",
      'cort' = "Corticosterone"
    )
               
    p <- ggplot(data = subset(d_sample, d_sample$treatment != ""), 
          mapping = aes(x = cap_doy)) + 
          geom_histogram(binwidth = 3.5, fill = "coral3", color = "gray30") + 
          theme_bw() + 
          facet_grid(treatment ~ group, labeller = labeller(treatment = as_labeller(treat_names), group = as_labeller(stage_names))) + 
          ylim(c(0, 13)) +
          ylab("Number of Samples") + xlab("Capture Day of Year") +
          geom_text(data = data.frame(treatment = c(rep("control", 3), rep("cort", 3)), 
                                      group = rep(c("1pre", "2post", "3cross"), 2),
                                      labels = LETTERS[1:6]), 
                    mapping = aes(x = -Inf, y = Inf, label = labels), hjust = -1, vjust = 1.5)
      
    ggsave(here("3_markdown_summary/samples.png"), p, device = "png", width = 8, height = 3.9)
                
# 0. Settings for MethylKit ----
    # Define the groups to be included
      group_1 <- "A"
      group_2 <- "D"
      
    # Define colors for plotting teh two groups (group 1 & group 2), see above in 'set colors' section for names
      col1 <- pre_col
      col2 <- cort_col
    
    # Define minimum coverage per site to be included
        min_cov <- 10
    
    # Define minimum SD for sites to be included for differential analysis
        min_sd <- 10   # removing sites that are essentially completely invariant
        min_sd_tile <- 10 
        
    # Define high percentage. Remove sites with percentage of reads at that site > than this. For PCR bias
        hi_perc <- 99.9
        
    # Define minimum number of sample per group for CpG to be included in comparison
        min_p_group <- 10L
        
    # Choose test or full (test runs only two samples as control treatment for speed)
        tester <- FALSE
        
    # Run with flat file database? Much slower but doesn't run out of memory
        db_yn <- TRUE
        
    # Use dss? (alternative is methylkit)
          dss <- "yes"
          if(dss == "no"){suff <- paste0(suff, "methylkit")}
          
    # Setting for tiles
          window_size <- 100
          step_size <- 100
          tile_cov <- 4       # Minimum cpgs in tile to include
          
    ## NO NEED TO CHANGE THESE
      # Labels for plots
          top_lab <- paste0(group_1, "vs", group_2)
      # Define suffix for naming any files or plots produced
          suff <- paste0("_", group_1, "v", group_2)
    
# 1. Read into methylkit ----
    # Make the file list for comparison
        # Define which groups to include in comparison (needs to be two comparisons)
            #d_sample$comp_grp <- d_sample$cross_grp
            sub_samples <- subset(d_sample, d_sample$comp_grp == group_1 | d_sample$comp_grp == group_2)
    
        # make a 1/0 coding for comparison groups
            ones <- data.frame(comp_grp = c(group_1, group_2), treat_code = c(0, 1))
            sub_samples <- plyr::join(sub_samples, ones, "comp_grp")
            
        if(tester == TRUE){
          sub_samples <- sub_samples[1:4, ]
          sub_samples$treat_code <- c(0, 0, 1, 1)
        }
    
    # Read into methylrawlist
          # No flat file database
            if(db_yn == FALSE){
              meth_list <- methRead(as.list(here("0_processed_data/bismark_cov_output", sub_samples$bis_processed_file)),
                                sample.id = as.list(sub_samples$sample_id),
                                assembly = "tres_20",
                                treatment = sub_samples$treat_code,
                                context = "CpG",
                                pipeline = "bismarkCoverage",
                                header = FALSE,
                                mincov = 10)
            }
            
          # Yes flat file database
            if(db_yn == TRUE){
              meth_list <- methRead(as.list(here("0_processed_data/bismark_cov_output", sub_samples$bis_processed_file)),
                                    sample.id = as.list(sub_samples$sample_id),
                                    assembly = "tres_20",
                                    treatment = sub_samples$treat_code,
                                    context = "CpG",
                                    pipeline = "bismarkCoverage",
                                    header = FALSE,
                                    dbtype = "tabix",
                                    dbdir = here("4_other_output"),
                                    mincov = 10)
            }
          
# 2. Filtering and uniting samples ----
        # Gives number of bases that fall into each bin of percent methylation
          # runs by sample. Could do as loop to make plots of all
            #getMethylationStats(meth_list[[5]], plot = TRUE)
            
        # Gives coverage per sample. Same as above.
            #getCoverageStats(meth_list[[5]], plot = TRUE)
            
    # Filters the list to exclude bases < x and > x coverage
        # greater than can be a problem because could be pcr bias
            f_meth_list <- filterByCoverage(meth_list, lo.count = min_cov, lo.perc = NULL,
                                            hi.count = NULL, hi.perc = hi_perc)
            
    # Combine the whole methyllist together into one object
        # normalize coverage
            f_meth_list2 <- normalizeCoverage(f_meth_list)
        # Set min per group to be minimum number of samples per group with coverage
            meth <- methylKit::unite(f_meth_list2, destrand = FALSE, min.per.group = min_p_group)

# 3. Save united object ----
    # Turn flat file database back into one object
      if(db_yn == TRUE){      
            meth2 <- as(meth, "methylBase")
      }
      if(db_yn == FALSE){
            meth2 <- meth
      }
    
    # Save the filtered, united, normalized object to RDS for reading in later                
        saveRDS(meth2, here("6_meth_RDS", paste0("meth", suff, ".RDS")))  
        meth2 <- readRDS(here("6_meth_RDS", paste0("meth", suff, ".RDS")))
        
# 3.5 START HERE IF READING IN ----        
      #meth2 <- readRDS(here("6_meth_RDS", "meth_AvD.rds"))   
            
# 4. Overall methylation percentage ----            
                      
    # make correlation matrix by sample of per CpG methylation
          pct_meth_cor <- cor(percMethylation(meth2), use = "complete.obs")
      
    # Plot correlation between samples. Only good if just a few samples.
          # getCorrelation(meth, plot = TRUE)
          
    # Vector of methylation percentages by site
          pct <- as.matrix(percMethylation(meth2))
          pc_cpg <- rep(NA, nrow(pct))
          pc_sd <- rep(NA, nrow(pct))
          pc_cv <- rep(NA, nrow(pct))
          for(i in 1:length(pc_cpg)){
            pc_cpg[i] <- mean(as.vector(pct[i, ]), na.rm = TRUE)
            pc_sd[i] <- sd(as.vector(pct[i, ]), na.rm = TRUE)
            pc_cv[i] <- pc_sd[i] / pc_cpg[i]
          }
        
      # Saves only CpG sites with SD for methylation higher than set value  
          include_list <- pc_sd > min_sd
          meth2 <- meth2[include_list, ]
          
# 5. Create PCA plot ----
        # Makes three panel PCA plot 1v2, 2v3, 1v3
        #plot clustered dendrogram
          clusterSamples(meth2, dist = "correlation", method = "ward", plot = TRUE)
          
        # PCA scree plot
          #PCASamples(meth, screeplot = TRUE)
          pc_meth <- PCASamples(meth2, obj.return = TRUE)
          
          pc_dat <- as.data.frame(pc_meth$x[, 1:3])
          pc_dat$sample_id <- rownames(pc_dat)
          pc_dat <- plyr::join(pc_dat, sub_samples, "sample_id", "left", "first")
          

          
          p1 <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC2, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC1") + ylab("Methylation PC2") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p2 <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC1") + ylab("Methylation PC3") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p3 <- ggplot(data = pc_dat, mapping = aes(x = PC2, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC2") + ylab("Methylation PC3") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            theme(legend.position = c(0.75, 0.11), legend.title = element_blank(), legend.background = element_blank())
          
          pa <- ggarrange(p1, p2, p3, nrow = 1)
          
          pa2 <- annotate_figure(pa,
                          top = text_grob(""),
                          fig.lab = top_lab,
                          fig.lab.pos = "top.left")
          ggsave(here("3_markdown_summary", paste0("pca", suff, ".png")), pa2, device = "png", width = 10, height = 3.6)
          
# 6. Finding differential CpGs ----
      # Use DSS beta-binomial model method to calculate differences       
          if(dss == "yes"){myDiff <- calculateDiffMethDSS(meth2, adjust = "SLIM")} 
          if(dss == "no"){myDiff <- calculateDiffMeth(meth2, test = "Chisq", overdispersion = "MN")}
        
      # save hyper methylated regions   (change to hypo for opposite)
          #diff25p_hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
          #diff25p_hypo <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hypo")
          
        # Plot percent methylation
          md <- as.data.frame(methylKit::getData(myDiff))
          md$pct_meth <- pc_cpg[include_list]
          p1 <- ggplot(data = md, mapping = aes(x = pc_cpg[include_list], y = -log10(qvalue))) + 
            theme_classic() +
            xlab("Overall Percent Methylation at CpG") +
            geom_smooth(col = "coral3", fill = "slateblue") +
            ylab("Average -log10(q-value)")
          
          cnt <- paste("Total number of CpGs", nrow(md), sep = " ")
          p2 <- ggplot(data = md, mapping = aes(pc_cpg[include_list])) + theme_classic() +
            geom_histogram(data = data.frame(pct = pc_cpg), mapping = aes(x = pct), fill = "slateblue", alpha = 0.5, binwidth = 2.5) +
            xlab("Overall Percent Methylation at CpG") +
            geom_histogram(fill = "coral3", color = "gray30", binwidth = 2.5, alpha = 0.5) +
            ylab("Number of CpG Sites") + ggtitle(cnt)
          
          p3 <- ggplot(data = md, mapping = aes(pc_cpg[include_list], y = abs(meth.diff))) + 
            theme_classic() +
            xlab("Overall Percent Methylation at CpG") +
            geom_smooth(col = "coral3", fill = "slateblue") +
            ylab("Abs(Group Difference in Methylation)")
          
          pb <- ggarrange(p1, p2, p3, nrow = 1)
          
          ggsave(here("3_markdown_summary", paste0("pct_diff", suff, ".png")), pb, device = "png", width = 10, height = 3.6)
            
         
        # Plot them, or turn off plot for text output 
           #dif_list <- diffMethPerChr(myDiff, plot = FALSE, qvalue.cutoff = 0.01, meth.cutoff = 25)
          
          p1 <- ggplot(data = methylKit::getData(myDiff), mapping = aes(x = pvalue)) +
            geom_histogram(fill = "coral3", binwidth = 0.01, color = "gray30") +
            theme_classic() + xlab("Difference p-value") +
            ylab("Number of CpG Sites")
          
          md2 <- md
          lt001 <- round(nrow(subset(md2, md2$qvalue < 0.05)))
          p2 <- ggplot(data = methylKit::getData(myDiff), mapping = aes(x = meth.diff, y = -log10(qvalue))) +
            geom_point(col = "slateblue", alpha = 0.7, size = 0.8) + theme_bw() +
            xlab("Difference in Methylation %") + ylab("-Log10(q-value)") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001)), linetype = "dashed") +
            ggtitle(paste(lt001, "CpGs q < 0.05", sep = " "))
          
          pc <- grid.arrange(p1, p2, nrow = 1)
          ggsave(here("3_markdown_summary", paste0("log_diff", suff, ".png")), pc, device = "png", width = 9, height = 4)
          
        # manhattan
          md2$seq <- seq(1, nrow(md2), 1)
          pm <- ggplot(data = md2, mapping = aes(x = seq, y = -log10(qvalue))) +
            theme_bw() + geom_point(alpha = 0.5, color = "slateblue", size = 0.9) + xlab("Position") + ylab("-Log10(q-value)") +
            geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001)), linetype = "dashed")
          ggsave(here("3_markdown_summary", paste0("manhattan", suff, ".png")), pm, device = "png", width = 10, height = 2.4)
          
  # Dendrogram
          #den <- clusterSamples(meth, dist = "correlation", method = "ward")
          #den$treatment <- sub_samples$trt_label
          #ggdendrogram(den)
          
  # Put plots together
          # ggsave(here("3_markdown_summary/xvx.png"),
          #   grid.arrange(pa2, pb, pc, nrow = 3),
          #   device = "png", width = 9.5, height = 10.2)
          
          
          
# 7. Create Methylkit Tiles ----
          # tile the genome to 1k bp windows and step across, object can be fed into unite and 
          # calculateDiffMeth functions consecutively to get differential regions
          # may want to use a lower coverage limit for this
          
          # read in again with lower minimum
          
          # tile
            meth_tile <- as(meth, "methylBase")
            tiles <- tileMethylCounts(meth_tile, win.size = window_size, step.size = step_size, cov.bases = tile_cov)
            tile2 <- normalizeCoverage(tiles)
            
            
          # Normalize tile coverage
            
          
          # if(db_yn == TRUE){      
          #   tiles2 <- tiles[,]
          # }
          # if(db_yn == FALSE){
          #   tiles2 <- tiles
          # }
          
# 8. TILE methylation percentage ----            
          
          # make correlation matrix by sample of per CpG methylation
          pct_meth_cor <- cor(percMethylation(tile2), use = "complete.obs")
          
          # Plot correlation between samples. Only good if just a few samples.
          # getCorrelation(meth, plot = TRUE)
          
          # Vector of methylation percentages by site
          pct <- as.matrix(percMethylation(tile2))
          pc_cpg <- rep(NA, nrow(pct))
          pc_sd <- rep(NA, nrow(pct))
          pc_cv <- rep(NA, nrow(pct))
          for(i in 1:length(pc_cpg)){
            pc_cpg[i] <- mean(as.vector(pct[i, ]), na.rm = TRUE)
            pc_sd[i] <- sd(as.vector(pct[i, ]), na.rm = TRUE)
            pc_cv[i] <- pc_sd[i] / pc_cpg[i]
          }
          
          # Saves only CpG sites with SD for methylation higher than set value  
          include_list <- pc_sd > min_sd_tile
          tile2 <- tile2[include_list, ]
          
# 5. TILE PCA plot ----
          # Makes three panel PCA plot 1v2, 2v3, 1v3
          #plot clustered dendrogram
          clusterSamples(tile2, dist = "correlation", method = "ward", plot = TRUE)
          
          # PCA scree plot
          #PCASamples(meth, screeplot = TRUE)
          pc_meth <- PCASamples(tile2, obj.return = TRUE)
          
          pc_dat <- as.data.frame(pc_meth$x[, 1:3])
          pc_dat$sample_id <- rownames(pc_dat)
          pc_dat <- plyr::join(pc_dat, sub_samples, "sample_id", "left", "first")
          
          
          
          p1 <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC2, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC1") + ylab("Methylation PC2") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p2 <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC1") + ylab("Methylation PC3") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p3 <- ggplot(data = pc_dat, mapping = aes(x = PC2, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC2") + ylab("Methylation PC3") +
            scale_fill_manual(values = c(col1, col2)) +
            scale_color_manual(values = c(col1, col2)) +
            theme(legend.position = c(0.75, 0.11), legend.title = element_blank(), legend.background = element_blank())
          
          pa <- ggarrange(p1, p2, p3, nrow = 1)
          
          pa2 <- annotate_figure(pa,
                                 top = text_grob(""),
                                 fig.lab = top_lab,
                                 fig.lab.pos = "top.left")
          ggsave(here("3_markdown_summary", paste0("TILE_pca", suff, ".png")), pa2, device = "png", width = 10, height = 3.6)
          
# 6. TILE differential CpGs ----
          # Use DSS beta-binomial model method to calculate differences       
          myDiff <- calculateDiffMethDSS(tile2, adjust = "SLIM")  
          
          # save hyper methylated regions   (change to hypo for opposite)
          #diff25p_hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
          #diff25p_hypo <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hypo")
          
          # Plot percent methylation
          md <- as.data.frame(getData(myDiff))
          md$pct_meth <- pc_cpg[include_list]
          p1 <- ggplot(data = md, mapping = aes(x = pc_cpg[include_list], y = -log10(qvalue))) + 
            theme_classic() +
            xlab("Overall Percent Methylation at CpG") +
            geom_smooth(col = "coral3", fill = "slateblue") +
            ylab("Average -log10(q-value)")
          
          cnt <- paste("Total number of CpGs", nrow(md), sep = " ")
          p2 <- ggplot(data = md, mapping = aes(pc_cpg[include_list])) + theme_classic() +
            geom_histogram(data = data.frame(pct = pc_cpg), mapping = aes(x = pct), fill = "slateblue", alpha = 0.5, binwidth = 2.5) +
            xlab("Overall Percent Methylation at CpG") +
            geom_histogram(fill = "coral3", color = "gray30", binwidth = 2.5, alpha = 0.5) +
            ylab("Number of CpG Sites") + ggtitle(cnt)
          
          p3 <- ggplot(data = md, mapping = aes(pc_cpg[include_list], y = abs(meth.diff))) + 
            theme_classic() +
            xlab("Overall Percent Methylation at CpG") +
            geom_smooth(col = "coral3", fill = "slateblue") +
            ylab("Abs(Group Difference in Methylation)")
          
          pb <- ggarrange(p1, p2, p3, nrow = 1)
          
          ggsave(here("3_markdown_summary", paste0("TILE_pct_diff", suff, ".png")), pb, device = "png", width = 10, height = 3.6)
          
          
          # Plot them, or turn off plot for text output 
          #dif_list <- diffMethPerChr(myDiff, plot = FALSE, qvalue.cutoff = 0.01, meth.cutoff = 25)
          
          p1 <- ggplot(data = getData(myDiff), mapping = aes(x = pvalue)) +
            geom_histogram(fill = "coral3", binwidth = 0.01, color = "gray30") +
            theme_classic() + xlab("Difference p-value") +
            ylab("Number of CpG Sites")
          
          md2 <- md
          lt001 <- round(nrow(subset(md2, md2$qvalue < 0.0001)) / nrow(md2) * 100, 2)
          p2 <- ggplot(data = getData(myDiff), mapping = aes(x = meth.diff, y = -log10(qvalue))) +
            geom_point(col = "slateblue", alpha = 0.7, size = 0.8) + theme_bw() +
            xlab("Difference in Methylation %") + ylab("-Log10(q-value)") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001)), linetype = "dashed") +
            ggtitle(paste(lt001, "% of Tiles q < 0.0001", sep = " "))
          
          pc <- grid.arrange(p1, p2, nrow = 1)
          ggsave(here("3_markdown_summary", paste0("TILE_log_diff", suff, ".png")), pc, device = "png", width = 9, height = 4)
          
          # manhattan
          md2$seq <- seq(1, nrow(md2), 1)
          pm <- ggplot(data = md2, mapping = aes(x = seq, y = -log10(qvalue))) +
            theme_bw() + geom_point(alpha = 0.5, color = "slateblue", size = 0.9) + xlab("Position") + ylab("-Log10(q-value)") +
            geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001)), linetype = "dashed")
          ggsave(here("3_markdown_summary", paste0("TILE_manhattan", suff, ".png")), pm, device = "png", width = 10, height = 2.4)
          
          # Dendrogram
          #den <- clusterSamples(meth, dist = "correlation", method = "ward")
          #den$treatment <- sub_samples$trt_label
          #ggdendrogram(den)
          
          # Put plots together
          # ggsave(here("3_markdown_summary/xvx.png"),
          #   grid.arrange(pa2, pb, pc, nrow = 3),
          #   device = "png", width = 9.5, height = 10.2)
          
          

# Other notes from MethylKit ----          
                    
    # Next section details how to put in covariates to the model to check for differential
          
    # Next section uses library(genomation) to see where DMRs are compared to annotated genome
          #e.g., what percent are in promoters, etc
          
    # Most of remainder of code is little tweaks and file type conversions, etc.
          
          # Section on batch effects to look at (and potentially remove) PCs associated with covariates
          #sampleAnnotation
          #assocComp
          #removeComp
          #percMethylation
          #reconstruct       
          
          
          
          
          
          
    
          
          
# DSS Approach ----
          
  # This is an alternative package that allows for using beta-binomial models.
      
# Pairwise correlation ----
    # Uses 'meth2' from above
          # Plot correlation between samples. Only good if just a few samples.
            cormat <- cor(percMethylation(meth2), use = "complete.obs")
            ltri <- lower.tri(cormat, diag = FALSE)
            cm2 <- reshape2::melt(cormat)
            lt2 <- reshape2::melt(ltri)
            colnames(lt2)[3] <- "include"
            cm2$include <- lt2$include
            cm2 <- subset(cm2, cm2$include == TRUE)
            
            ds_join <- d_sample[, c("sample_id", "band")]
            colnames(ds_join) <- c("Var1", "band1")            
            cm2 <- plyr::join(cm2, ds_join, "Var1", "left", "first")
            colnames(ds_join) <- c("Var2", "band2")
            cm2 <- plyr::join(cm2, ds_join, "Var2", "left", "first")
            for(i in 1:nrow(cm2)){
              ifelse(cm2$band1[i] == cm2$band2[i],
                     cm2$same[i] <- "Yes",
                     cm2$same[i] <- "No")
            }
            
            p_cor <- ggplot(data = cm2, mapping = aes(x = value, fill = same)) + geom_density() +
              scale_fill_viridis(discrete = TRUE, alpha = 0.5, begin = 0, end = 0.8) +
              theme_classic() + xlab("Correlation Coefficient") +
              ylab("Density") + guides(fill = guide_legend(title = "Repeat Sample")) + 
              theme(legend.position = c(0.2, 0.8)) + 
                guides(color = FALSE) +
                geom_rug(data = cm2, mapping = aes(color = same), sides = "b") +
                scale_color_viridis(discrete = TRUE, alpha = 0.5, begin = 0, end = 0.8)
            ggsave(here::here("3_markdown_summary", "pair_correlation.png"), p_cor, device = "png", width = 7.4, height = 6.4)
                             
     # 
              subs <- attr(meth2, which = "sample.ids")[1:5]
              sub_meth <- reorganize(meth2, sample.ids = subs, treatment = c(1, 1, 0, 0, 0))
              getCorrelation(sub_meth, plot = TRUE)
                             
# Meth sd by mean ----
     pct_sum <- data.frame(pc_cpg = pc_cpg, pc_sd = pc_sd, pc_cv = pc_cv) 
      pct_sum$Include <- "Yes"
      for(i in 1:nrow(pct_sum)){
        if(pct_sum$pc_sd[i] < 10){pct_sum$Include[i] <- "No"}
      }
    p_var <- ggplot(data = pct_sum, mapping = aes(x = pc_cpg, y = pc_sd, color = Include)) +
      geom_point(alpha = 0.1) + theme_classic() +
      scale_color_manual(values = c("slateblue", "coral3")) +
      geom_hline(yintercept = 10, linetype = "dashed") +
      xlab("Percent Methylation at CpG") +
      ylab("SD Between Samples") +
      theme(legend.position = c(0.1, 0.9)) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
    
    ggsave(here::here("3_markdown_summary", "pct_by_var.png"), p_var, device = "png", width = 7.4, height = 6.4)
    
## Methylation genome wide by sample ----
    pct <- as.matrix(percMethylation(meth2))
    pc_sample <- rep(NA, ncol(pct))
    for(i in 1:length(pc_sample)){
      pc_sample[i] <- mean(na.omit(pct[, i]))
    }
    pct_sam <- ggplot(data = data.frame(pc_sample = pc_sample), mapping = aes(x = pc_sample)) +
      geom_histogram(binwidth = 2, fill = "coral3", color = "gray30") +
      theme_classic() + xlab("Overall Percent Methylation") +
      ylab("Number of Samples")
    ggsave(here::here("3_markdown_summary", "pct_sample.png"), device = "png", width = 7.4, height = 6.4)      
              
    
# Arrows in PC plot ----
    pc_dat_pre <- subset(pc_dat, pc_dat$group_id == "1pre_cort")
    pc_dat_post <- subset(pc_dat, pc_dat$group_id == "2post_cort")
    post <- pc_dat_post[, c("band", "PC1", "PC2", "PC3")]
    colnames(post) <- c("band", "bPC1", "bPC2", "bPC3")    
    pc_dat_pre <- plyr::join(pc_dat_pre, post, "band")    
    pc_dat_pre$PC1_delta <- pc_dat_pre$PC1 - pc_dat_pre$bPC1
    pc_dat_pre$PC2_delta <- pc_dat_pre$PC2 - pc_dat_pre$bPC2
    pc_dat_pre$PC3_delta <- pc_dat_pre$PC3 - pc_dat_pre$bPC3
    
    p_delta <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC2, color = group_id)) + 
      geom_point(size = 4) + theme_bw() +
      scale_color_manual(values = c("gray70", "lightblue")) +
      guides(color = FALSE) +
      geom_segment(data = pc_dat_pre, mapping = aes(x = PC1, y = PC2, xend = bPC1, yend = bPC2),
                   arrow = arrow(type = "closed", length = unit(0.1, "in")))
    ggsave(here::here("3_markdown_summary", "delta_pc.png"), device = "png", width = 7.4, height = 6.4)

    
# Out to lme4 ----
    
    pct2 <- as.matrix(percMethylation(meth2))
    pc_cpg2 <- rep(NA, nrow(pct2))
    pc_sd2 <- rep(NA, nrow(pct2))
    pc_cv2 <- rep(NA, nrow(pct2))
    for(i in 1:length(pc_cpg2)){
      pc_cpg2[i] <- mean(as.vector(pct2[i, ]), na.rm = TRUE)
      pc_sd2[i] <- sd(as.vector(pct2[i, ]), na.rm = TRUE)
      pc_cv2[i] <- pc_sd2[i] / pc_cpg2[i]
    }
    
    
    pvals <- rep(NA, nrow(pct2))
    i <- 1
    for(i in 1:nrow(pct2)){
      sam <- colnames(pct2)
      pct_sam <- pct2[i, ]
      test <- data.frame(sample_id = sam, pct_sam = pct_sam)
      test2 <- plyr::join(test, pc_dat, "sample_id", "left", "first")
      m <- lmer(pct_sam ~ group_id + (1|band), data = test2)
      pvals[i] <- coefficients(summary(m))[2, 5]
    }
    
    hist(pvals, breaks = seq(0, 1, 0.01))
    
    saveRDS(pvals, here::here("5_temporary_files/pvals.rds"))
    pvals <- readRDS(here::here("5_temporary_files/pvals.rds"))
    
    plot(-log10(pvals), ylim = c(0, 6), xlim = c(0, 41000))
    abline(h = -log10(0.01), lty = 2, col = "coral3", lwd = 2)
    
    qobj <- qvalue::qvalue(p = pvals)
    
    hist(qobj)
    hist(qobj$qvalues)
    plot(-log10(qobj$qvalues), xlim = c(0, 41000))
    
  ## New approach for logistic regression  
    # Extracting info from methylkit (this could be a function)
      temp <- as.data.frame(methylKit::getData(meth2))
      temp2 <- temp %>%
        pivot_longer(cols = starts_with("coverage"),
                     names_to = "sample", 
                     values_to = "coverage")
      temp2 <- temp2[, c("chr", "start", "end", "strand", "sample", "coverage")]
      temp2$joiner <- paste(temp2$chr, temp2$start, temp2$end, temp2$sample, sep = "_")
      temp3 <- temp %>%
        pivot_longer(cols = starts_with("numCs"),
                     names_to = "sample2",
                     values_to = "num_cs")
      temp3 <- temp3[, c("chr", "start", "end", "sample2", "num_cs")]
      temp3$joiner <- paste(temp3$chr, temp3$start, temp3$end, temp3$sample2, sep = "_")
      temp3$joiner <- gsub("numCs", "coverage", temp3$joiner)
      temp4 <- plyr::join(temp2, temp3, "joiner", "left", "first")
      
      s_ids <- data.frame(sample_id = attr(meth2, which = "sample.ids"),
                          sample = paste0("coverage", seq(1, length(attr(meth2, which = "sample.ids")))))
      temp4 <- plyr::join(temp4, s_ids, "sample")
      
    # Joining to sample data
      ds <- d_sample[, c("sample_id", "band", "year", "min_age", "age_group", "treatment", 
                         "group", "group_id", "comp_grp", "trt_label", "cross_grp")]
    
    # Build data frame to store objects
      temp4$location <- paste(temp4$chr, temp4$start, temp4$end, sep = "_")

      
      #saveRDS(temp4, here::here("6_meth_RDS/extract_AvB.rds"))
      temp4 <- readRDS(here::here("6_meth_RDS/extract_AvB.rds"))
      t4 <- readRDS(here::here("6_meth_RDS/extract_DvE.rds"))
      
      u1 <- data.frame(location = unique(temp4$location))
      u2 <- data.frame(location = unique(t4$location), inAB = "yes")
      u1 <- plyr::join(u1, u2, "location", "left", "first")
      u3 <- subset(u1, u1$inAB == "yes")
      colnames(u3) <- c("locs", "inAB")
      
      temp4 <- temp4[, c("joiner", "chr", "start", "end", "sample", "coverage", "sample2", "num_cs", "sample_id", "location")]
      stemp4 <- filter(temp4, location %in% u3$locs)
      
      t4 <- t4[, c("joiner", "chr", "start", "end", "sample", "coverage", "sample2", "num_cs", "sample_id", "location")]
      st4 <- filter(t4, location %in% u3$locs)
      
      
      input <- rbind(st4, stemp4)
      output <- data.frame(location = unique(input$location),
                           n = NA, intercept = NA, grp_precort = NA, grp_postctl = NA, grp_postcort = NA,
                           reff_var = NA, dispersion = NA,
                           pAD = NA, pAB = NA, pAE = NA, pDB = NA, pDE = NA, pBE = NA, 
                           estAD = NA, estAB = NA, estAE = NA, estDB = NA, estDE = NA, estBE = NA,
                           seAD = NA, seAB = NA, seAE = NA, seDB = NA, seDE = NA, seBE = NA)
      output2 <- output
      
      #i <- 1
      
    # Loop through each cpg, run a model, and save output to the 'output' data frame
      for(i in 1:nrow(output)){
        # make a progress bar
          if(i == 1){pb <- txtProgressBar(min = 0, max = nrow(output), initial = 0, style = 3)}
        
        # make subset of data and join to metadata
          sub <- subset(input, input$location == output$location[i])    
          sub <- plyr::join(sub, ds, "sample_id")
        
        # fit model and extract emmeans contrasts
          m <- glmer(cbind(num_cs, coverage - num_cs) ~ group_id + (1|band), family = "binomial", data = sub,
                     control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8)))
          m2 <- glmmTMB(cbind(num_cs, coverage - num_cs) ~ group_id + (1|band), 
                                family = betabinomial(link = "logit"), 
                                data = sub)
          cps <- summary(contrast(emmeans(m, "group_id"), method = "pairwise", adjust = "none"))
          cps2 <- summary(contrast(emmeans(m2, "group_id"), method = "pairwise", adjust = "none"))
          
        # Calculate dispersion stat (in van oers code, from Zuur GLM & GLMM with R)
          output[i, "dispersion"] <- sum(residuals(m) ^ 2) / 
            (nrow(subset(sub, sub$coverage > 0)) - (length(fixef(m)) + 2))
          
          output2[i, "dispersion"] <- sum(residuals(m2) ^ 2) / 
            (nrow(subset(sub, sub$coverage > 0)) - (length(fixef(m2)) + 2))
          
        # Slot output into the saved data frame
          output[i, "n"] <- nrow(subset(sub, sub$coverage > 0))
          output[i, "intercept"] <- fixef(m)[1]
          output[i, "grp_precort"] <- fixef(m)[2]
          output[i, "grp_postctl"] <- fixef(m)[3]
          output[i, "grp_postcort"] <- fixef(m)[4]
          output[i, c("pAD", "pAB", "pAE", "pDB", "pDE", "pBE")] <- cps$p.value
          output[i, c("estAD", "estAB", "estAE", "estDB", "estDE", "estBE")] <- cps$estimate
          output[i, c("seAD", "seAB", "seAE", "seDB", "seDE", "seBE")] <- cps$SE
          output[i, "reff_var"] <- as.data.frame(VarCorr(m))[1, 5]
          
          output2[i, "n"] <- nrow(subset(sub, sub$coverage > 0))
          output2[i, "intercept"] <- fixef(m2)[[1]][1]
          output2[i, "grp_precort"] <- fixef(m2)[[1]][2]
          output2[i, "grp_postctl"] <- fixef(m2)[[1]][3]
          output2[i, "grp_postcort"] <- fixef(m2)[[1]][4]
          output2[i, c("pAD", "pAB", "pAE", "pDB", "pDE", "pBE")] <- cps2$p.value
          output2[i, c("estAD", "estAB", "estAE", "estDB", "estDE", "estBE")] <- cps2$estimate
          output2[i, c("seAD", "seAB", "seAE", "seDB", "seDE", "seBE")] <- cps2$SE
          #output2[i, "reff_var"] <- as.data.frame(VarCorr(m2))[1, 5]
        
        # update progress bar
          setTxtProgressBar(pb, i)
      }
      
      output$qAD <- qvalue::qvalue(output$pAD)$qvalues
      output$qAB <- qvalue::qvalue(output$pAB)$qvalues
      output$qDE <- qvalue::qvalue(output$pDE)$qvalues
      output$qBE <- qvalue::qvalue(output$pBE)$qvalues
      
      output2$qAD <- qvalue::qvalue(output2$pAD)$qvalues
      output2$qAB <- qvalue::qvalue(output2$pAB)$qvalues
      output2$qDE <- qvalue::qvalue(output2$pDE)$qvalues
      output2$qBE <- qvalue::qvalue(output2$pBE)$qvalues
      
      
      
      
      qs <- qvalue::qvalue(output$pvals)
      qs2 <- qvalue::qvalue(output$chisq)
      output$qvalue <- qs$qvalues
      output$qvalue_chi <- qs2$qvalues
      
      p1 <- ggplot(data = output2, mapping = aes(x = pAD)) + 
        geom_histogram(binwidth = 0.01, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "AvsD") + xlim(0, 1)
      p2 <- ggplot(data = output2, mapping = aes(x = pAB)) + 
        geom_histogram(binwidth = 0.01, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "AvsB") + xlim(0, 1)
      p3 <- ggplot(data = output2, mapping = aes(x = pDE)) + 
        geom_histogram(binwidth = 0.01, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "DvsE") + xlim(0, 1)
      p4 <- ggplot(data = output2, mapping = aes(x = pBE)) + 
        geom_histogram(binwidth = 0.01, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "BvsE") + xlim(0, 1)
      
      ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
      
      
      # Plot group differences
      difs <- subset(md, md$qvalue < 0.05)
      difs$loc <- paste(difs$chr, difs$start, difs$end, sep = "_")
      con <- data.frame(group_id = c("1pre_cort", "2post_cort"), grp2 = c(0, 1))
      #for(i in 1:39){
        i<-39
        ss <- subset(temp4, temp4$location == difs$loc[i])
        ss <- plyr::join(ss, ds, "sample_id")
        ss <- plyr::join(ss, con, "group_id", "left", "first")
        print(nrow(ss))
      #}
        
      dms <- filter(temp4, location %in% difs$loc)
      dms2 <- filter(t4, location %in% difs$loc)
      dms <- rbind(dms, dms2)
      dms <- plyr::join(dms, ds, "sample_id")
      
      dms$group_id <- factor(dms$group_id, levels = c("1pre_control", "2post_control", "1pre_cort", "2post_cort"))
    
    ggplot(data = ss, mapping = aes(x = group_id, y = num_cs / coverage)) + 
        geom_boxplot(fill = "orange", alpha = 0.6) + geom_point() + 
        geom_line(aes(group = as.factor(band)), color = "gray60", alpha = 0.5) +
        theme_classic() + xlab("") + ylab("Methylation Percentage")
    
    ggplot(data = dms, mapping = aes(x = group_id, y = num_cs / coverage, fill = group_id)) + 
      geom_point() +
      geom_boxplot(alpha = 0.6) +  
      scale_fill_manual(values = c(pre_col, cont_col, pre_col, cort_col)) +
      geom_line(aes(group = as.factor(band)), color = "gray60", alpha = 0.5) +
      theme_bw() + xlab("") + ylab("Methylation Percentage") +
      facet_wrap(~ location, scale = "free") +
      theme(strip.text = element_blank(), strip.background = element_blank(), axis.text.y = element_text(size = 6),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "bottom", legend.title = element_blank())
      
    p2 <- ggplot(data = ss, mapping = aes(x = grp2, y = num_cs / coverage, color = as.factor(band))) + 
      geom_point() + geom_line()
    ggarrange(p1, p2)
      
      