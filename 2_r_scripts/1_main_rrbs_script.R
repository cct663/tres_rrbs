# Script for tree swallow RRBS data from Ithaca NY
# Written by Conor Taff: cct663@gmail.com ~ cct63@cornell.edu
# Last updated 11/2020

# Largely based on vignette here: http://127.0.0.1:18146/library/methylKit/doc/methylKit.html

# Load libraries ----
    pacman::p_load("methylKit", "tidyverse", "here", "gridExtra", "ggpubr", "viridis", "DSS", "genomation",
                   "lme4", "lmerTest", "reshape2", "emmeans", "lmtest", "MuMIn", "qvalue")

# Load data ----
    cov_list <- list.files(here("0_processed_data/bismark_cov_output"))
    d_bis <- read.delim(here("0_processed_data/bismark_sample_summary.txt"))
    d_sample <- read.delim(here("1_raw_data/rrbs_sample_metadata.txt"))
    
# Set colors ----
    cort_col <- "#56B4E9"
    cont_col <- "#E69F00"
    pre_col <- "#999999"
    
    theme_rrbs <- function(){
      theme_bw() %+replace% # replace elements I want to change
        
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)
        ) 
    }
    
# Simple plots of reads and methylation ----
    # reads and methylation levels
            p_a <- d_bis %>%
                pivot_longer(cols = c("Total.Reads", "Aligned.Reads"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, breaks = seq(0, 20000, 1000), position = "identity", color = "gray30") +
                    theme_classic() + xlab("Sequences / 1000") +
                    scale_fill_manual(values = c(Total.Reads = "#F0E442", Aligned.Reads = "#56B4E9"),
                                      labels = c("Total Reads", "Aligned Reads")) +
                    ylab("Number of Samples") +
                    theme(legend.position = c(0.7, 0.8), legend.title = element_blank()) +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.7, vjust = 1.5, label = "A", size = 8) +
                    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 10), legend.text = element_text(size = 10))
            
            p_d <- d_bis %>%
                pivot_longer(cols = c("Methylated.CpGs", "Unmethylated.CpGs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, bins = 15, position = "identity", color = "gray30") +
                    theme_classic() + xlab("CpGs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CpGs = "#009E73", Methylated.CpGs = "#E69F00"),
                                      labels = c("Unmethylated", "Methylated")) +
                    ylab("Number of samples") +
                    #theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    guides(fill = "none") +
                    annotate("text", x = -Inf, y = Inf, hjust = -.9, vjust = 1.5, label = "D", size = 6) +
                    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 10))
            
            p_b <- d_bis %>%
                pivot_longer(cols = c("Methylated.CpHs", "Unmethylated.CpHs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, bins = 15, position = "identity", color = "gray30") +
                    theme_classic() + xlab("CpHs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CpHs = "#009E73", Methylated.CpHs = "#E69F00"),
                                      labels = c("Unmethylated", "Methylated")) +
                    ylab("Number of samples") + ylim(c(0, 130)) +
                    theme(legend.position = c(0.6, 0.65), legend.title = element_blank(), legend.text = element_text(size = 8)) +
                    annotate("text", x = -Inf, y = Inf, hjust = -1.5, vjust = 1.5, label = "B", size = 6) +
                    #annotate("text", x = 38e3, y = 25, label = "1.0% methylation", size = 3) +
                    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 10))
            
            p_c <- d_bis %>%
                pivot_longer(cols = c("Methylated.CHHs", "Unmethylated.CHHs"), names_to = "Type") %>%
                ggplot(mapping = aes(x = value/1000, fill = Type)) + 
                    geom_histogram(alpha = 0.5, bins = 15, position = "identity", color = "gray30") +
                    theme_classic() + xlab("CHHs / 1000") +
                    scale_fill_manual(values = c(Unmethylated.CHHs = "#009E73", Methylated.CHHs = "#E69F00"),
                                      labels = c("Unmethylated", "Methylated")) +
                    ylab("Number of samples") + ylim(c(0, 120)) +
                    #theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
                    guides(fill = "none") +
                    annotate("text", x = -Inf, y = Inf, hjust = -1.5, vjust = 1.5, label = "C", size = 6) +
                    #annotate("text", x = 75e3, y = 25, label = "0.6% methylation", size = 3) +
                    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 10)) +
                    scale_x_continuous(breaks = seq(0, 200000, 40000))
            
            p_e <- ggplot(data = d_bis, mapping = aes(x = (Methylated.CpGs / (Unmethylated.CpGs + Methylated.CpGs)) * 100)) +
                geom_histogram(fill = "#D55E00", binwidth = 2.5, color = "gray30", alpha = 0.6) + 
                theme_classic() + xlab("CpG methylation %") +
                ylab("Number of samples") +
                annotate("text", x = -Inf, y = Inf, hjust = -.9, vjust = 1.5, label = "E", size = 6) +
                theme(axis.text = element_text(size = 11), axis.title = element_text(size = 10))
            
            ggsave(here("3_markdown_summary/read_summary.png"),
                grid.arrange(p_a, p_b, p_c, p_d, p_e,
                         layout_matrix = rbind(c(1, 1, 2, 4),
                                               c(1, 1, 3, 5))),
                device = "png", width = 9, height = 4.5)
    
    
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
      group_1 <- "C"
      group_2 <- "F"
      
    # Define colors for plotting teh two groups (group 1 & group 2), see above in 'set colors' section for names
      col1 <- pre_col
      col2 <- cort_col
    
    # Define minimum coverage per site to be included
        min_cov <- 10
    
    # Define minimum SD for sites to be included for differential analysis
        min_sd <- 5   # removing sites that are essentially completely invariant
        
    # Define high percentage. Remove sites with percentage of reads at that site > than this. For PCR bias
        hi_perc <- 99.5
        
    # Define minimum number of sample per group for CpG to be included in comparison
        min_p_group <- 6L  # changing to 6 for cross year, at 10 for within year
        
    # Randomization number of loops (takes a while)
        n_rand <- 10
        
    # Choose test or full (test runs only two samples as control treatment for speed)
        tester <- FALSE
        
    # Run with flat file database? Much slower but doesn't run out of memory
        db_yn <- TRUE
        
    # Use dss? (alternative is methylkit)
          dss <- "yes"
          if(dss == "no"){suff <- paste0(suff, "methylkit")}
          
    ## NO NEED TO CHANGE THESE
      # Labels for plots
          top_lab <- paste0(group_1, "vs", group_2)
      # Define suffix for naming any files or plots produced
          suff <- paste0("_", group_1, "v", group_2)
    
# 1. Read into methylkit ----
    # Make the file list for comparison
        # Define which groups to include in comparison (needs to be two comparisons)
            #d_sample$comp_grp <- d_sample$cross_grp  # use for cross year only
            sub_samples <- subset(d_sample, d_sample$comp_grp == group_1 | d_sample$comp_grp == group_2)
    
        # make a 1/0 coding for comparison groups
            ones <- data.frame(comp_grp = c(group_1, group_2), treat_code = c(0, 1))
            sub_samples <- plyr::join(sub_samples, ones, "comp_grp")
            
        if(tester == TRUE){
          sub_samples <- sub_samples[1:4, ]
          sub_samples$treat_code <- c(0, 0, 1, 1)
        }
            
        # change to character
            sub_samples$sample_id <- as.character(sub_samples$sample_id)
    
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
            
          # Yes flat file database. Fits more in memory but takes longer.
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
        #saveRDS(meth2, here("6_meth_RDS", paste0("meth", suff, ".RDS")))  
        #meth2 <- readRDS(here("6_meth_RDS", paste0("meth", suff, ".RDS")))
        
# 3.5 START HERE IF READING IN ----        
      #meth2 <- readRDS(here::here("6_meth_RDS", "meth_AvD.rds")) 
          # for basic description of patterns use meth_AvD
            
# 4. Overall methylation percentage ----            
                      
    # make correlation matrix by sample of per CpG methylation
          pct_meth_cor <- cor(percMethylation(meth2), use = "complete.obs")
      
    # Plot correlation between samples. Only good if just a few samples.
          # getCorrelation(meth, plot = TRUE)
          
    # Vector of methylation percentages by site
          pct <- as.matrix(percMethylation(meth2))
          pc_cpg <- rep(NA, nrow(pct))
          pc_median <- rep(NA, nrow(pct))
          pc_sd <- rep(NA, nrow(pct))
          pc_cv <- rep(NA, nrow(pct))
          for(i in 1:length(pc_cpg)){
            pc_cpg[i] <- mean(as.vector(pct[i, ]), na.rm = TRUE)
            pc_median[i] <- median(as.vector(pct[i, ]), na.rm = TRUE)
            pc_sd[i] <- sd(as.vector(pct[i, ]), na.rm = TRUE)
            pc_cv[i] <- pc_sd[i] / pc_cpg[i]
          }
        
      # Saves only CpG sites with SD for methylation higher than set value  
          il1 <- pc_sd < min_sd   # gets rid of very low variation sites
          il2 <- pc_median == 0 | pc_median == 100   # gets rid of basically invariant sites
          ild <- data.frame(il1 = il1, il2 = il2)
          ild$il3 <- ild$il1 + ild$il2
          include_list <- ild$il3 == 0
          meth2 <- meth2[include_list, ]
          
          # pct <- as.matrix(percMethylation(meth2))  # running again after filtering low sd
          # pct <- as.data.frame(pct)     # for test at end run with filtered meth2
          # pct$chr <- meth2$chr
          # pct$start <- meth2$start
          # pct$end <- meth2$end
          # pct$strand <- meth2$strand
          # pct_ad <- pct
          
          
          
# 5. Create PCA plot ----
        ## NOT REPORTED
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
            #scale_fill_manual(values = c(col1, col2)) +
            #scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p2 <- ggplot(data = pc_dat, mapping = aes(x = PC1, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC1") + ylab("Methylation PC3") +
            #scale_fill_manual(values = c(col1, col2)) +
            #scale_color_manual(values = c(col1, col2)) +
            guides(fill = FALSE, color = FALSE) +
            theme(legend.position = c(0.85, 0.9), legend.title = element_blank(), legend.background = element_blank())
          
          p3 <- ggplot(data = pc_dat, mapping = aes(x = PC2, y = PC3, fill = group_id, color = group_id)) +
            geom_point(pch = 21, size = 3, alpha = 0.7) + 
            theme_bw() + stat_ellipse() + 
            xlab("Methylation PC2") + ylab("Methylation PC3") +
            #scale_fill_manual(values = c(col1, col2)) +
            #scale_color_manual(values = c(col1, col2)) +
            theme(legend.position = c(0.75, 0.11), legend.title = element_blank(), legend.background = element_blank())
          
          pa <- ggarrange(p1, p2, p3, nrow = 1)
          
          pa2 <- annotate_figure(pa,
                          top = text_grob(""),
                          fig.lab = top_lab,
                          fig.lab.pos = "top.left")
          ggsave(here("3_markdown_summary", paste0("pca", suff, ".png")), pa2, device = "png", width = 10, height = 3.6)
          
# 6. Finding differential CpGs ----
      ## THIS IS BASIC METHYLKIT APPROACH BUT NOT USING BECAUSE DOESNT FIT OUR DATA STRUCTURE    
          
      # Use DSS beta-binomial model method to calculate differences       
          if(dss == "yes"){myDiff <- calculateDiffMethDSS(meth2, adjust = "fdr", mc.cores = 4)} 
          if(dss == "no"){myDiff <- calculateDiffMeth(meth2, test = "F", adjust = "SLIM", overdispersion = "MN", mc.cores = 4)}
        
      # save hyper methylated regions   (change to hypo for opposite)
          #diff25p_hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
          #diff25p_hypo <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hypo")
          diff_all <- getMethylDiff(myDiff, qvalue = 0.05, type = "all", difference = 0)
          
        # Plot percent methylation
          md <- as.data.frame(methylKit::getData(myDiff))
          md$pct_meth <- pc_cpg[include_list]
          p1 <- ggplot(data = md, mapping = aes(x = pc_cpg[include_list], y = -log10(qvalue))) + 
            theme_classic() +
            xlab("Overall Percent Methylation at CpG") +
            geom_smooth(col = "coral3", fill = "slateblue") +
            ylab("Average -log10(q-value)")
          
          cnt <- paste("Total CpGs", nrow(md), sep = " ")
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
          
          
          plot_data <- methylKit::getData(myDiff)
          pd0t10 <- subset(plot_data, plot_data$meth.diff > -10 & plot_data$meth.diff < 10)
          pd10p <- subset(plot_data, plot_data$meth.diff > 10 | plot_data$meth.diff < -10)
          p1 <- ggplot(data = plot_data, mapping = aes(x = pvalue)) +
            geom_histogram(fill = "coral3", binwidth = 0.01, color = "gray30", boundary = 0) +
            theme_classic() + xlab("Difference p-value") +
            ylab("Number of CpG Sites")
          
          md2 <- md
          lt001 <- round(nrow(subset(md2, md2$qvalue < 0.05)))
          p2 <- ggplot(data = plot_data, mapping = aes(x = meth.diff, y = -log10(pvalue))) +
            geom_point(col = "slateblue", alpha = 0.7, size = 0.8) + theme_bw() +
            xlab("Difference in Methylation %") + ylab("-Log10(q-value)") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001)), linetype = "dashed") +
            ggtitle(paste(lt001, "CpGs q < 0.05", sep = " "))
          
          pc <- grid.arrange(p1, p2, nrow = 1)
          ggsave(here("3_markdown_summary", paste0("log_diff", suff, ".png")), pc, device = "png", width = 9, height = 4)
          
        # manhattan
          md2$seq <- seq(1, nrow(md2), 1)
          pm <- ggplot(data = md2, mapping = aes(x = seq, y = -log10(pvalue))) +
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
          
# 7. Randomization ----
          ## NOT USING THIS CURRENTLY 
      # 
      # i <- 1
      # trts <- attr(meth2, which = "treatment") 
      # meth_r <- meth2
      # r_save <- data.frame(loop = seq(1, n_rand, 1), r_05 = NA, r_01 = NA, r_001 = NA)
      # 
      # for(i in 1:n_rand){
      #     # make a progress bar
      #       if(i == 1){pb <- txtProgressBar(min = 0, max = n_rand, initial = 0, style = 3)}
      #     trts_r <- trts[order(runif(length(trts), 0, 1))]    
      #     
      #     
      #     attr(meth_r, which = "treatment") <- trts_r
      #     
      #     if(dss == "yes"){r_diff <- calculateDiffMethDSS(meth_r, adjust = "SLIM", mc.cores = 4)} 
      #     if(dss == "no"){r_diff <- calculateDiffMeth(meth_r, test = "Chisq", adjust = "hochberg", overdispersion = "MN", mc.cores = 4)}
      #     
      #     md_r <- as.data.frame(methylKit::getData(r_diff))
      #     r_save[i, "r_05"] <- round(nrow(subset(md_r, md_r$pvalue < 0.05)))
      #     r_save[i, "r_01"] <- round(nrow(subset(md_r, md_r$pvalue < 0.01)))
      #     r_save[i, "r_001"] <- round(nrow(subset(md_r, md_r$pvalue < 0.001)))
      #     
      #     # update progress bar
      #       setTxtProgressBar(pb, i)
      # }    
      # 
      # 
      #   p1 <- ggplot(data = r_save, mapping = aes(x = r_05)) + 
      #     geom_histogram(fill = "gray70", color = "gray30", binwidth = 2, alpha = 0.6) +
      #     xlab("CpGs with q < 0.05") + ylab("Count") +
      #     theme_classic() +
      #     geom_vline(xintercept = round(nrow(subset(md2, md2$pvalue < 0.05))), linetype = "dashed", color = "coral3", size = 2)
      #   
      #   p2 <- ggplot(data = r_save, mapping = aes(x = r_01)) + 
      #     geom_histogram(fill = "gray70", color = "gray30", binwidth = 2, alpha = 0.6) +
      #     xlab("CpGs with q < 0.01") + ylab("Count") +
      #     theme_classic() +
      #     geom_vline(xintercept = round(nrow(subset(md2, md2$pvalue < 0.01))), linetype = "dashed", color = "coral3", size = 2)
      #   
      #   p3 <- ggplot(data = r_save, mapping = aes(x = r_001)) + 
      #     geom_histogram(fill = "gray70", color = "gray30", binwidth = 2, alpha = 0.6) +
      #     xlab("CpGs with q < 0.001") + ylab("Count") +
      #     theme_classic() +
      #     geom_vline(xintercept = round(nrow(subset(md2, md2$pvalue < 0.001))), linetype = "dashed", color = "coral3", size = 2)
      #       
      #    
      # 
      # 
      # pr <- grid.arrange(p1, p2, p3, nrow = 1)
      # ggsave(here("3_markdown_summary", paste0("randomization", suff, ".png")), pc, device = "png", width = 9.5, height = 3.5)
      #     

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

    
# Out to GLMM models ----
    
  ## Approach for logistic regression of each CPG site 
    # Extracting info from methylkit (this could be a function). Need to run separately for pairs then combine
      # Three comparisons cover everything: AvD (pre), BvE (post), FvC (year later) 
      # Run three times and save to file with right suffix (using code chunks above)
        # 2/10/23 this is already run and saved to file to load in below
              # temp <- as.data.frame(methylKit::getData(meth2))  # meth2 should be one of the comparisons AD/BE/FC
              # temp2 <- temp %>%
              #   pivot_longer(cols = starts_with("coverage"),
              #                names_to = "sample",
              #                values_to = "coverage")
              # temp2 <- temp2[, c("chr", "start", "end", "strand", "sample", "coverage")]
              # temp2$joiner <- paste(temp2$chr, temp2$start, temp2$end, temp2$sample, sep = "_")
              # temp3 <- temp %>%
              #   pivot_longer(cols = starts_with("numCs"),
              #                names_to = "sample2",
              #                values_to = "num_cs")
              # temp3 <- temp3[, c("chr", "start", "end", "sample2", "num_cs")]
              # temp3$joiner <- paste(temp3$chr, temp3$start, temp3$end, temp3$sample2, sep = "_")
              # temp3$joiner <- gsub("numCs", "coverage", temp3$joiner)
              # temp4 <- plyr::join(temp2, temp3, "joiner", "left", "first")
              # 
              # s_ids <- data.frame(sample_id = attr(meth2, which = "sample.ids"),
              #                     sample = paste0("coverage", seq(1, length(attr(meth2, which = "sample.ids")))))
              # temp4 <- plyr::join(temp4, s_ids, "sample")
              # temp4 <- na.omit(temp4)

        # Joining to sample data
          # ds <- d_sample[, c("sample_id", "band", "year", "min_age", "age_group", "treatment",
          #                    "group", "group_id", "comp_grp", "trt_label", "cross_grp",
          #                    "mass", "bhead", "fwing", "s_cort", "b_cort")]
          # temp4 <- plyr::join(temp4, ds, "sample_id")


      # save the produced object. saves for each of the three are loaded from file
          #saveRDS(temp4, here::here("6_meth_RDS/extract_FvC.rds"))   # change suffix to match comparison
          filt_DA <- readRDS(here::here("6_meth_RDS/extract_DvA.rds")) #DA is 1st capture year 1
          filt_BE <- readRDS(here::here("6_meth_RDS/extract_BvE.rds")) #BE is post treatment capture year 1
          filt_CF <- readRDS(here::here("6_meth_RDS/extract_FvC.rds")) #CF is capture in 2016 one year after treatments
      
      # add in location column
          filt_DA$location <- paste(filt_DA$chr, filt_DA$start, filt_DA$end, sep = "_")
          filt_BE$location <- paste(filt_BE$chr, filt_BE$start, filt_BE$end, sep = "_")
          filt_CF$location <- paste(filt_CF$chr, filt_CF$start, filt_CF$end, sep = "_")
      
      # filter to locations in both datasets for pre-post within year
        # this is just figuring out which CpGs actually have data for both pre and post
            u1 <- data.frame(location = unique(filt_DA$location))   # list of cpg locations in D&A group
            u2 <- data.frame(location = unique(filt_BE$location), inDA = "yes") # list of cpg locations in B&E group
            u1 <- plyr::join(u1, u2, "location", "left", "first") #joins list of sites in DA vs. BE together
            u3 <- subset(u1, u1$inDA == "yes") # removes sites from list that aren't in both objects
            colnames(u3) <- c("locs", "inDA")
            
            filt_DA <- filt_DA[, c("joiner", "chr", "start", "end", "sample", "coverage", 
                                   "sample2", "num_cs", "sample_id", "location",
                                   "band", "year", "min_age", "age_group", "treatment", 
                                   "group", "comp_grp", "trt_label", "mass", "bhead", "fwing", "s_cort", "b_cort")]
            stemp4 <- filter(filt_DA, location %in% u3$locs)
            
            filt_BE <- filt_BE[, c("joiner", "chr", "start", "end", "sample", "coverage", 
                                   "sample2", "num_cs", "sample_id", "location",
                                   "band", "year", "min_age", "age_group", "treatment", 
                                   "group", "comp_grp", "trt_label", "mass", "bhead", "fwing", "s_cort", "b_cort")]
            st4 <- filter(filt_BE, location %in% u3$locs)
            
      # filter to locations in both datasets for cross year comparison. same as last code block but for cross year
          # this is just figuring out which CpGs actually have data for yr 1 and yr 2
            u1x <- data.frame(location = unique(filt_DA$location))
            u2x <- data.frame(location = unique(filt_CF$location), inDA = "yes")
            u1x <- plyr::join(u1x, u2x, "location", "left", "first")
            u3x <- subset(u1x, u1x$inDA == "yes")
            colnames(u3x) <- c("locs", "inDA")
            
            filt_DA <- filt_DA[, c("joiner", "chr", "start", "end", "sample", "coverage", 
                                   "sample2", "num_cs", "sample_id", "location",
                                   "band", "year", "min_age", "age_group", "treatment", 
                                   "group", "comp_grp", "trt_label", "mass", "bhead", "fwing", "s_cort", "b_cort")]
            stemp4x <- filter(filt_DA, location %in% u3x$locs)
            
            filt_CF <- filt_CF[, c("joiner", "chr", "start", "end", "sample", "coverage", 
                                   "sample2", "num_cs", "sample_id", "location",
                                   "band", "year", "min_age", "age_group", "treatment", 
                                   "group", "comp_grp", "trt_label", "mass", "bhead", "fwing", "s_cort", "b_cort")]
            st4x <- filter(filt_CF, location %in% u3x$locs)
        
      # filter for the correlation with natural cort models
          # in this case don't need the same joining as done above
            nat_cort <- filt_DA[, c("joiner", "chr", "start", "end", "sample", "coverage", 
                                    "sample2", "num_cs", "sample_id", "location",
                                    "band", "year", "min_age", "age_group", "treatment", 
                                    "group", "comp_grp", "trt_label", "mass", "bhead", "fwing", "s_cort", "b_cort")]
            
            
          # Add column for percent methylation  
            stemp4$pct_meth <- stemp4$num_cs / stemp4$coverage  # within year pre
            st4$pct_meth <- st4$num_cs / st4$coverage # within year post
            nat_cort$pct_meth <- nat_cort$num_cs / nat_cort$coverage # unmanipulated cort comparison
            
            stemp4x$pct_meth <- stemp4x$num_cs / stemp4x$coverage # cross year pre
            st4x$pct_meth <- st4x$num_cs / st4x$coverage # cross year post
            
            stemp4$pp_join <- paste(stemp4$location, stemp4$band, sep = "_") # within year pre
            st4$pp_join <- paste(st4$location, st4$band, sep = "_") # within year post
            
            stemp4x$pp_join <- paste(stemp4x$location, stemp4x$band, sep = "_") # cross year pre
            st4x$pp_join <- paste(st4x$location, st4x$band, sep = "_") # cross year post
            
          # join reduced version of pre measures to the post measures
            # first for within year comparison
              stemp4b <- stemp4[, c("pp_join", "coverage", "num_cs", "pct_meth", "sample_id",
                              "group", "mass", "bhead", "fwing", "s_cort", "b_cort")]
              colnames(stemp4b) <- c("pp_join", "pre_cov", "pre_numcs", "pre_pctm", "pre_sampid",
                                  "pre_group", "pre_mass", "pre_bhead", "pre_fwing", "pre_scort", "pre_bcort")
              comb_pp <- plyr::join(st4, stemp4b, "pp_join")
              
            # next for cross year comparison
              stemp4bx <- stemp4x[, c("pp_join", "coverage", "num_cs", "pct_meth", "sample_id",
                                    "group", "mass", "bhead", "fwing", "s_cort", "b_cort")]
              colnames(stemp4bx) <- c("pp_join", "pre_cov", "pre_numcs", "pre_pctm", "pre_sampid",
                                     "pre_group", "pre_mass", "pre_bhead", "pre_fwing", "pre_scort", "pre_bcort")
              comb_ppx <- plyr::join(st4x, stemp4bx, "pp_join")

      # make wider version for pre-post comparison
          # gets rid of post rows where no pre measure is available for within season
            comb_pp <- na.omit(comb_pp) # not sure of effect on sample sizes for comps?loses 25% reads
            comb_ppb <- comb_pp[, c("coverage", "num_cs", "location", "band", "treatment",
                                    "pre_pctm")]  
            
          # gets rid of post rows where no pre measure is available for cross season
            comb_ppx <- na.omit(comb_ppx) # not sure of effect on sample sizes for comps?loses 25% reads
            comb_ppbx <- comb_ppx[, c("coverage", "num_cs", "location", "band", "treatment",
                                    "pre_pctm")] 
      
      # prepare an output file to store model results from loops
          # for pre-post
            output <- data.frame(location = unique(comb_ppb$location),
                               n = NA, trt_p = NA, trt_eff = NA,
                               pre_p = NA, pre_eff = NA, intercept = NA,
                               singular = NA, r2m = NA, r2c = NA, message = NA)
            
          # for cross year
            outputx <- data.frame(location = unique(comb_ppbx$location),
                                  n = NA, trt_p = NA, trt_eff = NA,
                                  pre_p = NA, pre_eff = NA, intercept = NA,
                                  singular = NA, r2m = NA, r2c = NA, message = NA)
            
          # for initial capture correlation of cort and methylation
            output_cor <- data.frame(location = unique(nat_cort$location),
                                     base_n = NA, base_int = NA, base_eff = NA, base_p = NA,
                                     base_singular = NA, base_r2c = NA, base_r2m = NA, base_message = NA,
                                     si_n = NA, si_int = NA, si_eff = NA, si_p = NA,
                                     si_singular = NA, si_r2c = NA, si_r2m = NA, si_message = NA)  
              
      # fit models for each cpg
            # first for within year comparison  
                start <- Sys.time()
                  #note i = 156 fails useful for testing
                  for(i in 1:nrow(output)){
                    sub <- subset(comb_ppb, comb_ppb$location == output$location[i])
                    m <- glmer(cbind(num_cs, coverage - num_cs) ~ treatment + scale(pre_pctm) + (1|band),
                                           family = "binomial", data = sub,
                                            control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2,
                                                                  optCtrl = list(maxfun = 2e8)))
                       
                    if(length(m@optinfo$conv$lme4$messages) > 0){
                      output$message[i] <- m@optinfo$conv$lme4$messages
                    } 
                    
                    output$n[i] <- nrow(sub)
                    
                    
                    s <- summary(m)
                    output$trt_p[i] <- s$coefficients[2, 4]
                    output$trt_eff[i] <- s$coefficients[2, 1]
                    output$pre_p[i] <- s$coefficients[3, 4]
                    output$pre_eff[i] <- s$coefficients[3, 1]
                    output$intercept[i] <- s$coefficients[1,1]
                    
                    output$singular[i] <- isSingular(m)
                    
                    rs <- suppressWarnings(r.squaredGLMM(m))
                    output$r2m[i] <- rs[1,1]
                    output$r2c[i] <- rs[1,2]
                    
                    em_m <- pairs(emmeans(m, ~ treatment), adjust = "none")
                    
                    print(i)
                  }
                end <- Sys.time()
                    within_elapse <- end - start
                    
                    output$cont_est <- rethinking::logistic(output$intercept)
                    output$cort_est <- rethinking::logistic(output$intercept + output$trt_eff)
                    output$cont_min_cort <- output$cont_est - output$cort_est
                    qq <- qvalue(output$trt_p, fdr.level = 0.05)  
                    output$trt_q <- qq$qvalues
                    qq_pre <- qvalue(output$pre_p, fdr.level = 0.05)
                    output$pre_q <- qq_pre$qvalues
                    output2 <- subset(output, is.na(output$message) == TRUE)
                    output2$sig <- "no"
                    for(i in 1:nrow(output2)){
                      if(output2$trt_q[i] < 0.05){output2$sig[i] <- "yes"}
                    }
                    
                    sig_list <- subset(output2, output2$trt_q < 0.05)
          
            # Next for cross year comparison
                    start <- Sys.time()
                    for(i in 1:nrow(outputx)){
                      subx <- subset(comb_ppbx, comb_ppbx$location == outputx$location[i])
                      outputx$n[i] <- nrow(subx)
                      if(outputx$n[i] > 6){
                          mx <- glmer(cbind(num_cs, coverage - num_cs) ~ treatment + scale(pre_pctm) + (1|band),
                                     family = "binomial", data = subx,
                                     control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2,
                                                            optCtrl = list(maxfun = 2e8)))
                          
                          if(length(mx@optinfo$conv$lme4$messages) > 0){
                            outputx$message[i] <- mx@optinfo$conv$lme4$messages
                          } 
                          
                          
                          
                          
                          sx <- summary(mx)
                          outputx$trt_p[i] <- sx$coefficients[2, 4]
                          outputx$trt_eff[i] <- sx$coefficients[2, 1]
                          outputx$pre_p[i] <- sx$coefficients[3, 4]
                          outputx$pre_eff[i] <- sx$coefficients[3, 1]
                          outputx$intercept[i] <- sx$coefficients[1,1]
                          
                          outputx$singular[i] <- isSingular(mx)
                          
                          rsx <- suppressWarnings(r.squaredGLMM(mx))
                          outputx$r2m[i] <- rsx[1,1]
                          outputx$r2c[i] <- rsx[1,2]
                          
                          em_mx <- pairs(emmeans(mx, ~ treatment), adjust = "none")
                      }
                      
                      print(i)
                    }
                    end <- Sys.time()
                    cross_elapse <- end - start
                    
                    outputx$cont_est <- rethinking::logistic(outputx$intercept)
                    outputx$cort_est <- rethinking::logistic(outputx$intercept + outputx$trt_eff)
                    outputx$cont_min_cort <- outputx$cont_est - outputx$cort_est
                    qqx <- qvalue(outputx$trt_p, fdr.level = 0.05)  
                    qq_prex <- qvalue(outputx$pre_p, fdr.level = 0.05)
                    outputx$trt_q <- qqx$qvalues
                    outputx$pre_q <- qq_prex$qvalues
                    
                    outputx$sig <- "no"
                    for(i in 1:nrow(outputx)){
                      if(is.na(outputx$trt_q[i]) == FALSE){
                        if(outputx$trt_q[i] < 0.05){outputx$sig[i] <- "yes"}
                      }
                    }
                    output2x <- subset(outputx, outputx$singular == FALSE)
                    output2x <- subset(output2x, is.na(output2x$message) == TRUE)
                    
                    
                    sig_listx <- subset(output2x, output2x$trt_q < 0.05)
                    
            # Finally for pre-treatment correlation with natural cort variation
                  
                    start <- Sys.time()
                    for(i in 1:nrow(output_cor)){
                      subc <- subset(nat_cort, nat_cort$location == output_cor$location[i])
                      output_cor$base_n[i] <- nrow(subset(subc, is.na(subc$b_cort) == FALSE))
                      output_cor$si_n[i] <- nrow(subset(subc, is.na(subc$s_cort) == FALSE))
                      
                      if(output_cor$base_n[i] > 9){
                        mc <- glmer(cbind(num_cs, coverage - num_cs) ~ scale(b_cort) + (1|band),
                                  family = "binomial", data = subc)
                        
                        if(length(mc@optinfo$conv$lme4$messages) > 0){
                          output_cor$base_message[i] <- mc@optinfo$conv$lme4$messages
                        } 

                        sc <- summary(mc)$coefficients
                        rc <- suppressWarnings(r.squaredGLMM(mc))
                        output_cor$base_r2c[i] <- rc[1, 2]
                        output_cor$base_r2m[i] <- rc[1, 1]
                        output_cor$base_int[i] <- sc[1, 1]
                        output_cor$base_eff[i] <- sc[2, 1]
                        output_cor$base_p[i] <- sc[2, 4]
                        output_cor$base_singular[i] <- isSingular(mc)
                      }
                      
                      if(output_cor$si_n[i] > 9){
                        ms <- glmer(cbind(num_cs, coverage - num_cs) ~ scale(s_cort) + (1|band),
                                    family = "binomial", data = subc)
                        
                        if(length(ms@optinfo$conv$lme4$messages) > 0){
                          output_cor$si_message[i] <- ms@optinfo$conv$lme4$messages
                        } 
                        
                        ss <- summary(ms)$coefficients
                        rs <- suppressWarnings(r.squaredGLMM(ms))
                        output_cor$si_r2c[i] <- rs[1, 2]
                        output_cor$si_r2m[i] <- rs[1, 1]
                        output_cor$si_int[i] <- ss[1, 1]
                        output_cor$si_eff[i] <- ss[2, 1]
                        output_cor$si_p[i] <- ss[2, 4]
                        output_cor$si_singular[i] <- isSingular(ms)
                      }
                    
                      
                      
                      
                      print(i)
                    }
                    end <- Sys.time()
                    obs_elapse <- end - start 
                    
                    qqb <- qvalue(output_cor$base_p, fdr.level = 0.05)  
                    output_cor$base_q <- qqb$qvalues
                    qqs <- qvalue(output_cor$si_p, fdr.level = 0.05)
                    output_cor$si_q <- qqs$qvalues
                    
                    output_cor$base_sig <- "no"
                    output_cor$si_sig <- "no"
                    
                    for(i in 1:nrow(output_cor)){
                      if(output_cor$base_q[i] < 0.05){output_cor$base_sig[i] <- "yes"}
                      if(output_cor$si_q[i] < 0.05){output_cor$si_sig[i] <- "yes"}
                    }
                    
                    output_base <- subset(output_cor, output_cor$base_singular == FALSE & is.na(output_cor$base_message) == TRUE)
                    output_si <- subset(output_cor, output_cor$si_singular == FALSE & is.na(output_cor$si_message) == TRUE)
                    
                    
                    sig_list_base <- subset(output_base, output_base$base_q < 0.05)
                    sig_list_si <- subset(output_si, output_si$si_q < 0.05)
                    
          # collecting the output of all of these loops into one object to save
              # so that it can be loaded without running again
                    
                model_loop_output <- list(output, output2, sig_list,     # within year
                                            outputx, output2x, sig_listx,  # between year
                                            output_cor,                    # correlation with cort
                                            output_base, sig_list_base,    # base cort only
                                            output_si, sig_list_si)       # induced cort only
                saveRDS(model_loop_output, here::here("4_other_output/model_loop_output.rds"))
                model_loop_output <- readRDS(here::here("4_other_output/model_loop_output.rds"))
          
          # make a volcano plot for the four comparisons
              output_base <- model_loop_output[[8]]
              pv1 <- ggplot(data = output_base, 
                            mapping = aes(x = base_eff, y = -log(as.numeric(base_p), 10), color = base_sig)) +
                geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001), 10), linetype = "dashed", color = "gray70") +
                theme_rrbs() +
                geom_point(alpha = 0.7, size = 0.8) +
                xlab("Regression coefficient of \n baseline corticosterone") +
                ylab("-log10(p-value)") +
                scale_color_manual(values = c("gray20", "red")) +
                guides(color = "none") +
                coord_cartesian(xlim = c(-3, 3), ylim = c(0, 6.5)) +
                ggtitle("Baseline corticosterone") +
                #ggtitle("116 of 78,143 CpGs") + # associated with baseline corticosterone (78,143)") +
                annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
                theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
              
              output_si <- model_loop_output[[10]]
              pv2 <- ggplot(data = output_si, 
                            mapping = aes(x = si_eff, y = -log(as.numeric(si_p), 10), color = si_sig)) +
                geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001), 10), linetype = "dashed", color = "gray70") +
                theme_rrbs() +
                geom_point(alpha = 0.7, size = 0.8) +
                xlab("Regression coefficient of \n stress-induced corticosterone") +
                ylab("-log10(p-value)") +
                scale_color_manual(values = c("gray20", "red")) +
                guides(color = "none") +
                coord_cartesian(xlim = c(-3, 3), ylim = c(0, 6.5)) +
                ggtitle("Stress-induced corticosterone") +
                #ggtitle("356 of 78,027 CpGs") + # associated with stress-induced corticosterone (78,027)") +
                annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
                theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
              
              within_output <- model_loop_output[[2]]
              pv3 <- ggplot(data = within_output, 
                            mapping = aes(x = cont_min_cort * 100, y = -log(as.numeric(trt_p), 10), color = sig)) +
                geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001), 10), linetype = "dashed", color = "gray70") +
                theme_rrbs() +
                geom_point(alpha = 0.7, size = 0.8) +
                xlab("Difference in % methylation \n control minus treatment") +
                ylab("-log10(p-value)") +
                scale_color_manual(values = c("gray20", "red")) +
                coord_cartesian(ylim = c(0, 6.5), xlim = c(-100, 100)) +
                guides(color = "none") +
                ggtitle("Within year treatment") +
                #ggtitle("111 of 48,070 CpGs") + # associated with treatment (48,070)") +
                annotate(geom = "text", label = "C", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
                theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
              
              across_output <- model_loop_output[[5]]
              pv4 <- ggplot(data = across_output, 
                            mapping = aes(x = cont_min_cort * 100, y = -log(as.numeric(trt_p), 10), color = sig)) +
                geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001), 10), linetype = "dashed", color = "gray70") +
                theme_rrbs() +
                geom_point(alpha = 0.7, size = 0.8) +
                xlab("Difference in % methylation \n control minus treatment") +
                ylab("-log10(p-value)") +
                scale_color_manual(values = c("gray20", "red")) +
                coord_cartesian(ylim = c(0, 6.5), xlim = c(-100, 100)) +
                guides(color = "none") +
                ggtitle("Between year treatment") +
                #ggtitle("49 of 6,787 CpGs") + # associated with treatment (6,787)") +
                annotate(geom = "text", label = "D", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
                theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
              
              
              pviol <- ggpubr::ggarrange(pv1, pv2, pv3, pv4, ncol = 2, nrow = 2)
              saveRDS(pviol, file = here::here("5_temporary_files/pviol.rds"))
              
              #ggsave(here::here("2_r_scripts/pviol.png"), plot = pviol,
              #       device = "png", width = 7.8, height = 6.6, units = "in", dpi = 300)
              
      # make a similar plot for pre treatment methylation percentage
          pre_w <- ggplot(data = output2, mapping = aes(x = pre_eff, y = -log(pre_p, 10))) +
            theme_rrbs() +
            xlab("Regression coefficient \n of pre- vs. post-treatment \n methylation %") +
            ylab("-log10(p-value)") +
            annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
            geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001)), linetype = "dashed", color = "gray70") +
            geom_point(size = 0.8, alpha = 0.7, color = "gray20") +
            ggtitle("Within year sampling") +
            coord_cartesian(xlim = c(-5.5, 5.5), ylim = c(0, 32)) +
            theme(axis.title.x = element_text(size = 13))
          
          pre_b <- ggplot(data = output2x, mapping = aes(x = pre_eff, y = -log(pre_p, 10))) +
            theme_rrbs() +
            xlab("Regression coefficient \n of pre- vs. post-treatment \n methylation %") +
            ylab("-log10(p-value)") +
            annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
            geom_hline(yintercept = -log(c(0.05, 0.01, 0.001, 0.0001)), linetype = "dashed", color = "gray70") +
            geom_point(size = 0.8, alpha = 0.7, color = "gray20") +
            ggtitle("Between year sampling") +
            coord_cartesian(xlim = c(-5.5, 5.5), ylim = c(0, 32)) +
            theme(axis.title.x = element_text(size = 13))
          
          consistent_meth <- ggarrange(pre_w, pre_b, nrow = 1)
          saveRDS(consistent_meth, here::here("5_temporary_files/consistent_meth.rds"))
            
              
      # use the significant lists to filter down methylkit object so I can map to genomic features
            # start by reading in methylkit for AvD at line way above as 'meth2'
              # get out the location names from meth2 and see which ones match significant list
                methx <- methylKit::getData(meth2)[, 1:7]
                methx$location <- paste(methx$chr, methx$start, methx$end, sep = "_")
                # this one is for within year treatment
                  sl_within <- sig_list[, c("location", "sig")]
                  methx2 <- data.frame(location = methx$location, rownum = seq(1, nrow(methx), 1))
                  sl_within <- plyr::join(sl_within, methx2, "location")
                  within_sig <- meth2[sl_within$rownum, ]
                
                # this one is for between year treatment
                  sl_between <- sig_listx[, c("location", "sig")]
                  methx2b <- data.frame(location = methx$location, rownum = seq(1, nrow(methx), 1))
                  sl_between <- plyr::join(sl_between, methx2b, "location")
                  between_sig <- meth2[sl_between$rownum, ]
                
                # this one is for baseline corticosterone
                  sl_base <- sig_list_base[, c("location", "base_sig")]
                  methx2c <- data.frame(location = methx$location, rownum = seq(1, nrow(methx), 1))
                  sl_base <- plyr::join(sl_base, methx2c, "location")
                  base_sig <- meth2[sl_base$rownum, ]  
                  
                # this one is for induced corticosterone
                  sl_si <- sig_list_si[, c("location", "si_sig")]
                  methx2d <- data.frame(location = methx$location, rownum = seq(1, nrow(methx), 1))
                  sl_si <- plyr::join(sl_si, methx2d, "location")
                  si_sig <- meth2[sl_si$rownum, ]
                
                  
                             
          
          
### BELIEVE THIS IS ALL OLD BELOW HERE IN THIS SECTION      ----
                  
      # for cort
      #       output_c <- data.frame(location = unique(input$location),
      #                              n = NA, intercept = NA, b_cort_beta = NA, s_cort_beta = NA, dispersion = NA,
      #                              b_cort_lrt = NA, s_cort_lrt = NA,
      #                              sd_meth = NA, mu_meth = NA)
      #     
      # # for cross season
      #       input_x <- subset(t4_CAD, is.na(t4_CAD$num_cs) == FALSE)
      #       input_x$pct_meth <- input_x$num_cs / input_x$coverage
      #       output_x <- data.frame(location = unique(input_x$location),
      #                              n = NA, grp_cntl = NA, grp_cort = NA, dispersion = NA, x_lrt = NA,
      #                              sd_meth = NA, mu_meth = NA)

      
      # Find sites that don't differ at all between any groups (not enough change)
            # for within season
              input_s <- plyr::join(input, ds, "sample_id")
              difs <- input_s %>%
                group_by(location, group_id) %>%
                summarise(pc_meth = mean(pct_meth)) %>%
                pivot_wider(id_cols = location, names_from = group_id, values_from = pc_meth, names_prefix = "w")
              difs$maxdif <- NA
              for(i in 1:nrow(difs)){
                difs$maxdif[i] <- max(difs[i, 2:5]) - min(difs[i, 2:5])
              }
              difs$ppc_dif <- abs(difs$w2post_cort + difs$w2post_control - difs$w1pre_cort - difs$w1pre_control)
              difs2 <- as.data.frame(difs[, c("location", "maxdif", "ppc_dif")])
              output <- plyr::join(output, difs2, "location")
              output_c <- plyr::join(output_c, difs2, "location")
              
            # for between season
              input_xs <- plyr::join(input_x, ds, "sample_id")
              difs_x <- input_xs %>%
                group_by(location, cross_grp) %>%
                summarise(pc_meth = mean(pct_meth)) %>%
                pivot_wider(id_cols = location, names_from = cross_grp, values_from = pc_meth, names_prefix = "w")
              difs_x$diff <- abs(difs_x$wF - difs_x$wCAD)
              difs_x2 <- as.data.frame(difs_x[, c("location", "diff")])
              output_x <- plyr::join(output_x, difs_x2, "location")
      
      #mall <- lmer(pct_meth ~ 1 + (1|location) + (1|sample_id), data = input)
    
    # Loop for pre-post and for cort 
    # Loop through each cpg, run a model, and save output to the 'output' data frame
      start <- Sys.time()
      for(i in 1:nrow(output)){
        # make a progress bar
          if(i == 1){pb <- txtProgressBar(min = 0, max = nrow(output), initial = 0, style = 3)}
        
        # make subset of data and join to metadata
          sub <- subset(input, input$location == output$location[i])    
          sub <- plyr::join(sub, ds, "sample_id")
          sub2 <- subset(sub, sub$group == "1pre" & is.na(sub$b_cort) == FALSE & is.na(sub$s_cort) == FALSE) # used just for cort correlations
          
        # # figure out correlation between repeat samples from the same individual
        #   subw <- sub %>%
        #     pivot_wider(id_cols = band, names_from = group, names_prefix = "w", values_from = pct_meth)
        #   subw <- na.omit(subw)
        
        # fit model
          m <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ treatment*group + (1|band) + (1|sample_id), family = "binomial", data = sub,
                     control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
          m_b <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ treatment + group + (1|band) + (1|sample_id), family = "binomial", data = sub,
                                        control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
          output[i, "inter_lrt"] <- lrtest(m, m_b)[2, 5]
          output[i, "sd_meth"] <- sd(sub$pct_meth)
          output[i, "mu_meth"] <- mean(sub$pct_meth)
          output$n[i] <- nrow(sub)
          output$grp_prectl[i] <- exp(fixef(m)[1])
          output$grp_precort[i] <- exp(fixef(m)[1] + fixef(m)[2])
          output$grp_postctl[i] <- exp(fixef(m)[1] + fixef(m)[3])
          output$grp_postcort[i] <- exp(fixef(m)[1] + fixef(m)[2] + fixef(m)[4])
          
          em_m <- as.data.frame(pairs(emmeans(m, ~ treatment*group), adjust = "none"))
          output$comp_AB_e[i] <- em_m[2, 2]
          output$comp_AD_e[i] <- em_m[1, 2]
          output$comp_AE_e[i] <- em_m[3, 2]
          output$comp_BD_e[i] <- em_m[4, 2]
          output$comp_BE_e[i] <- em_m[6, 2]
          output$comp_DE_e[i] <- em_m[5, 2]
          
          output$comp_AB_p[i] <- em_m[2, 6]
          output$comp_AD_p[i] <- em_m[1, 6]
          output$comp_AE_p[i] <- em_m[3, 6]
          output$comp_BD_p[i] <- em_m[4, 6]
          output$comp_BE_p[i] <- em_m[6, 6]
          output$comp_DE_p[i] <- em_m[5, 6]
          
        # Calculate dispersion stat (in van oers code, from Zuur GLM & GLMM with R)
          output[i, "dispersion"] <- sum(residuals(m) ^ 2) /
            (nrow(sub) - (length(fixef(m)) + 2))
          
        # fit model for cort pre treatment
          mc <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ scale(b_cort) + scale(s_cort) + (1|sample_id), 
                                        family = "binomial", data = sub2, control = glmerControl(optimizer = "bobyqa", 
                                        boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
          mc1 <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ scale(s_cort) + (1|sample_id), 
                                       family = "binomial", data = sub2, control = glmerControl(optimizer = "bobyqa", 
                                        boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
          mc2 <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ scale(b_cort) + (1|sample_id), 
                                       family = "binomial", data = sub2, control = glmerControl(optimizer = "bobyqa", 
                                        boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
        
        # fill in output for cort
            output_c[i, "sd_meth"] <- sd(sub2$pct_meth)
            output_c[i, "mu_meth"] <- mean(sub2$pct_meth)
            output_c[i, "n"] <- nrow(sub2)
            output_c[i, "intercept"] <- fixef(mc)[1]
            output_c[i, "b_cort_beta"] <- fixef(mc)[2]
            output_c[i, "s_cort_beta"] <- fixef(mc)[3]
            output_c[i, "b_cort_lrt"] <- lrtest(mc, mc1)[2, 5]
            output_c[i, "s_cort_lrt"] <- lrtest(mc, mc2)[2, 5]
            output_c[i, "dispersion"] <- sum(residuals(mc) ^ 2) /
              (nrow(sub2) - (length(fixef(mc)) + 2))
          
        
        # update progress bar
          setTxtProgressBar(pb, i)
          #print(i)
      }
      end <- Sys.time()
      end - start
      
    # Loop and add in correlation in pre to post samples for each cpg
      output$pp_cor <- NA
      output_c$pp_cor <- NA
      output$repeatability <- NA
      for(i in 1:nrow(output)){
        
        # make a progress bar
        #if(i == 1){pb <- txtProgressBar(min = 0, max = nrow(output_x), initial = 0, style = 3)}
        
        # make subset of data and join to metadata
          sub <- subset(input, input$location == output$location[i])    
          sub <- plyr::join(sub, ds, "sample_id")
          
          an <- anova(lm(pct_meth ~ as.factor(band), data = sub))
          output$repeatability[i] <- ((an$`Mean Sq`[1] - an$`Mean Sq`[2]) / 2) / 
            (((an$`Mean Sq`[1] - an$`Mean Sq`[2]) / 2) + an$`Mean Sq`[2])
        
        # figure out correlation between repeat samples from the same individual
            # subw <- sub %>%
            #   pivot_wider(id_cols = band, names_from = group, names_prefix = "w", values_from = pct_meth)
            # subw <- na.omit(subw)
        
        # add correlation
            #output$pp_cor[i] <- cor(subw$w1pre, subw$w2post)
            
            # update progress bar
            #setTxtProgressBar(pb, i)
            print(i)
      }

      
      # loop for cross year
        for(i in 1:nrow(output_x)){
          # make a progress bar
            if(i == 1){pb <- txtProgressBar(min = 0, max = nrow(output_x), initial = 0, style = 3)}
          
          # make subset of data and join to metadata
            sub <- subset(input_x, input_x$location == output_x$location[i])    
            sub <- plyr::join(sub, ds, "sample_id")
               
          # fit model and extract emmeans contrasts for pre post comparison
            mx <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ cross_grp + (1|band) + (1|sample_id), family = "binomial", data = sub,
                                        control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
            m_bx <- suppressMessages(glmer(cbind(num_cs, coverage - num_cs) ~ 1 + (1|band) + (1|sample_id), family = "binomial", data = sub,
                                          control = glmerControl(optimizer = "bobyqa", boundary.tol = 1e-2, optCtrl = list(maxfun = 2e8))))
            output_x[i, "x_lrt"] <- lrtest(mx, m_bx)[2, 5]
            output_x[i, "sd_meth"] <- sd(sub$pct_meth)
            output_x[i, "mu_meth"] <- mean(sub$pct_meth)
            output_x$n[i] <- nrow(sub)
            output_x$grp_cntl[i] <- exp(fixef(mx)[1])
            output_x$grp_cort[i] <- exp(fixef(mx)[1] + fixef(mx)[2])
          
          # Calculate dispersion stat (in van oers code, from Zuur GLM & GLMM with R)
            output_x[i, "dispersion"] <- sum(residuals(mx) ^ 2) /
              (nrow(sub) - (length(fixef(mx)) + 2))
          
          # update progress bar
            setTxtProgressBar(pb, i)
            #print(i)
        }
      
      end <- Sys.time()
      end - start
      
      # for pre-post
        #saveRDS(output, here::here("4_other_output/output.RDS"))
        #output <- readRDS(here::here("4_other_output/output.RDS"))
      
        outtest <- output %>%
          filter(location %in% outs$location)
        
      # for cort
        #saveRDS(output_c, here::here("4_other_output/output_c.RDS"))
        #output_c <- readRDS(here::here("4_other_output/output_c.RDS"))
        
      # for cort
        #saveRDS(output_x, here::here("4_other_output/output_x.RDS"))
        #output_x <- readRDS(here::here("4_other_output/output_x.RDS"))
      
      output$q_inter_lrt <- qvalue::qvalue(output$inter_lrt)$qvalues
      output_c$q_st_p <- qvalue::qvalue(output_c$s_cort_lrt)$qvalues
      
      #output <- plyr::join(output, difs2, "location")
      outputd <- subset(output, output$repeatability < 0.5 &
                          output$mu_meth > 0.1 & output$mu_meth < 0.9 &
                          output$ppc_dif > 0.1)
      hist(outputd$inter_lrt, breaks = seq(0, 1, 0.01))
      nrow(outputd)
      outputd$pa_inter <- p.adjust(outputd$inter_lrt, method = "fdr")
      nrow(subset(outputd, outputd$pa_inter < 0.05))
      hist(outputd$pa_inter, breaks = seq(0, 1, 0.01))
      outputd$q_inter_lrt <- qvalue::qvalue(outputd$inter_lrt)$qvalues
      nrow(subset(outputd, outputd$q_inter_lrt < 0.05))
      hist(outputd$q_inter_lrt, breaks = seq(0, 1, 0.01))
      oo <- qvalue::qvalue(outputd$inter_lrt)
      
      outputd_x <- subset(output_x, output_x$mu_meth > 0.2 & output_x$mu_meth < 0.8 &
                            output_x$sd_meth/output_x$mu_meth > 0.4)
      outputd_x$pa_inter <- p.adjust(outputd_x$x_lrt, method = "fdr")
      outputd_x$q_x_lrt <- qvalue::qvalue(outputd_x$x_lrt)$qvalues
      nrow(outputd_x)
      hist(outputd_x$x_lrt, breaks = seq(0, 1, 0.01))
      hist(outputd_x$q_x_lrt, breaks = seq(0, 1, 0.01))
      nrow(subset(outputd_x, outputd_x$q_x_lrt < 0.05))
      hist(outputd_x$pa_inter, breaks = seq(0, 1, 0.01))
      nrow(subset(outputd_x, outputd_x$pa_inter < 0.05))
      
      outputd_x <- subset(output_x, output_x$diff > 0.1)
      
      
      # wide
        in2 <- pivot_wider(input, id_cols = sample_id, names_from = location, values_from = pct_meth)
        nb <- estim_ncpPCA(in2[, 2:ncol(in2)], method.cv = "Kfold", verbose = FALSE)
        in2na <- imputePCA(in2[, 2:ncol(in2)], ncp = nb$ncp)
        in3 <- in2[, colSums(is.na(in2)) == 0]
        prin <- prcomp(in3[2:ncol(in3)], center = TRUE)
        prin2 <- data.frame(sample_id = in3$sample_id, pr1 = prin$x[, 1], pr2 = prin$x[, 2], pr3 = prin$x[, 3])
        prin2 <- plyr::join(prin2, ds, "sample_id")
        ggplot(prin2, mapping = aes(x = pr1, y = pr2, color = group_id)) + 
          geom_point(size = 2) +
          theme_bw() +
          stat_ellipse()
        mxx <- lmer(pr1 ~ treatment*group + (1|band), data = prin2)
      
      # To plot groups from a subset
          sub2 <- pivot_wider(sub, id_cols = c(band, treatment), names_from = group, values_from = pct_meth, 
                              names_prefix = "t_")
          ggplot(sub2, mapping = aes(x = t_1pre, y = t_2post, color = treatment)) +
            geom_smooth(method = "lm") +
            geom_jitter()
        
          sub$xplot <- as.numeric(sub$group) - 1 
          ggplot(data = sub, mapping = aes(x = xplot, y = pct_meth, color = treatment, by = as.factor(band))) +
            geom_point() +
            geom_line() +
            geom_boxplot(mapping = aes(by = as.factor(xplot), y = pct_meth, color = treatment), inherit.aes = FALSE)
          
          ggplot(data = sub, mapping = aes(x = group, y = pct_meth)) +
            geom_violin() +
            geom_jitter(width = 0.2)
            
      
      qs <- qvalue::qvalue(output$pvals)
      qs2 <- qvalue::qvalue(output$chisq)
      output$qvalue <- qs$qvalues
      output$qvalue_chi <- qs2$qvalues
      
      p1 <- ggplot(data = output2, mapping = aes(x = pAD)) + 
        geom_histogram(binwidth = 0.05, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "AvsD") + xlim(0, 1)
      p2 <- ggplot(data = output2, mapping = aes(x = pAB)) + 
        geom_histogram(binwidth = 0.05, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "AvsB") + xlim(0, 1)
      p3 <- ggplot(data = output2, mapping = aes(x = pDE)) + 
        geom_histogram(binwidth = 0.05, fill = "coral3", color = "gray30") +
        theme_classic() + xlab("q-value") + ggtitle(label = "DvsE") + xlim(0, 1)
      p4 <- ggplot(data = output2, mapping = aes(x = pBE)) + 
        geom_histogram(binwidth = 0.05, fill = "coral3", color = "gray30") +
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
      
      
    

    
# Genomic features ----
    gene_ob <- readTranscriptFeatures(here::here("7_annotated_genome/tres_annotated.bed"),
                                                  remove.unusual = FALSE,
                                                  up.flank = 2000,
                                                  down.flank = 0)
    
    # Annotate differences to genomic features
      gene_an <- suppressWarnings(annotateWithGeneParts(as(meth2, "GRanges"), gene_ob)) # all pre-treatment
      #gene_an2 <- annotate.WithGenicParts(methx_a2, gene_ob)
      
    # find percent methylation by region type
      # meth2_sub <- meth2[gene_an@dist.to.TSS$target.row, ]
      # ppct <- as.matrix(percMethylation(meth2_sub))
      # ppc_cpg <- rep(NA, nrow(ppct))
      # for(i in 1:length(ppc_cpg)){
      #   ppc_cpg[i] <- mean(as.vector(ppct[i, ]), na.rm = TRUE)
      # }
      # region_pct <- data.frame(pct = ppc_cpg,
      #                          row = gene_an@dist.to.TSS$target.row,
      #                          region = "intergenic")
      # check <- gene_an@members
      # for(i in 1:nrow(region_pct)){
      #   if(check[i, 3] == 1){region_pct$region[i] <- "intron"}
      #   if(check[i, 2] == 1){region_pct$region[i] <- "exon"}
      #   if(check[i, 1] == 1){region_pct$region[i] <- "promoter"}
      # }
      # 
      # ggplot(region_pct, mapping = aes(x = pct, y = region)) +
      #   geom_density(adjust = 5)
      
    # get nearest gene name
      gene_match <- getAssociationWithTSS(gene_an)
      
      
    # regional analysis
      promoters <- regionCounts(meth2, gene_ob$promoters)
      exons <- regionCounts(meth2, gene_ob$exons)
      introns <- regionCounts(meth2, gene_ob$introns)
      tssers <- regionCounts(meth2, gene_ob$TSSes)
      
    # make a df of annotated locations
      df_annot <- data.frame(chr = c(as.character(promoters@.Data[[1]]), as.character(exons@.Data[[1]]), as.character(introns@.Data[[1]]), as.character(tssers@.Data[[1]])),
                             start = c(promoters@.Data[[2]], exons@.Data[[2]], introns@.Data[[2]], tssers@.Data[[2]]),
                             end = c(promoters@.Data[[3]], exons@.Data[[3]], introns@.Data[[3]], tssers@.Data[[3]]),
                             type = c(rep("promoter", nrow(promoters)), rep("exon", nrow(exons)), rep("intron", nrow(introns)), rep("tss", nrow(tssers))))
      df_annot$location <- paste(df_annot$chr, df_annot$start, df_annot$end, sep = "_")
      
    # join
        hh <- plyr::join(df_annot, output, "location", type = "left")
        oox <- output
        for(i in 1:nrow(oox)){
          sps <- str_split(oox$location[i], "_")
          oox$chr[i] <- sps[[1]][1]
          oox$start[i] <- as.numeric(sps[[1]][2])
          oox$end[i] <- as.numeric(sps[[1]][3])
        }
        oox$type <- NA
        for(i in 1:nrow(output_x)){
          sub <- subset(df_annot, df_annot$chr == oox$chr[i] & df_annot$start <= oox$start[i] & df_annot$end >= oox$end[i])
          if(nrow(sub) > 1){
            oox$type[i] <- as.character(sub$type[1])
          }
        }
        ggplot(oox, mapping = aes(x = inter_lrt)) + 
          geom_histogram(binwidth = 0.01, boundary = 0, fill = "coral3") +
          facet_wrap(~ type, scales = "free")
        
        oox2 <- subset(oox, oox$type == "promoter")
    
    # Pct meth for regions
      sd(na.omit(percMethylation(promoters)))
      sd(na.omit(percMethylation(exons)))
      sd(na.omit(percMethylation(introns)))
      
    # summarize to make plot of methylation by genomic feature
      pmeth <- rowMeans(as.matrix(percMethylation(promoters)), na.rm = TRUE)
      emeth <- rowMeans(as.matrix(percMethylation(exons)), na.rm = TRUE)
      imeth <- rowMeans(as.matrix(percMethylation(introns)), na.rm = TRUE)
      
      df_reg <- data.frame(pct = c(pmeth, emeth, imeth),
                           reg = c(rep("promoter", length(pmeth)),
                                   rep("exon", length(emeth)),
                                   rep("intron", length(imeth))))
      
      dfreg2 <- df_reg %>%
        group_by(reg) %>%
        summarise(mu = mean(pct), med = median(pct), sem = sd(pct) / sqrt(n()))
      sum2 <- ggplot(data = df_reg, mapping = aes(x = reg, y = pct)) +
        geom_boxplot(notch = TRUE, fill = "lightblue") +
        theme_rrbs() +
        xlab("") +
        ylab("Percent methylation") +
        geom_point(data = dfreg2, inherit.aes = FALSE, mapping = aes(x = reg, y = mu)) +
        annotate(geom = "text", label = "B", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
        theme(axis.text.x = element_text(size = 11))
        #theme(axis.text.x = element_text(angle = 30, hjust = 0.95))
      
      sum1 <- ggplot(as.data.frame(pc_cpg), mapping = aes(x = pc_cpg)) +
        geom_histogram(breaks = seq(0, 100, 4), fill = "lightblue", color = "gray40") +
        theme_rrbs() +
        coord_cartesian(ylim = c(0, 52000)) +
        annotate(geom = "text", label = "A", x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5, size = 7) +
        ylab("Number of CpGs") +
        xlab("Percent methylation")
      
      sum_plot <- ggpubr::ggarrange(sum1, sum2, widths = c(2, 1))
      saveRDS(sum_plot, here::here("5_temporary_files/sum_plot.rds"))
      
      # ggplot(data = df_reg, mapping = aes(x = pct)) +
      #   geom_histogram(bindwidth = 4) +
      #   facet_wrap(~ reg, ncol = 1)
      
      
      
    # Association with tss
      diffAnn <- getAssociationWithTSS(gene_an)
      targAnn <- getTargetAnnotationStats(gene_an, percentage = TRUE, precedence = TRUE)
      plotTargetAnnotation(gene_an, precedence = TRUE, main = "differential")
      getFeatsWithTargetsStats(gene_an, percentage = TRUE)
      
## testing removing based on similarity in repeat samples DELETE BELOW HERE?? ----
    
      #taking from line 303 output
      
      pct_ad2 <- pct_ad %>%
        as.data.frame() %>% 
        pivot_longer(cols = T_045:T_089, names_to = "sample_id") %>%
        na.omit()
      
      pct_be2 <- pct_be %>%
        as.data.frame() %>%
        pivot_longer(cols = RR_006:RR_036, names_to = "sample_id") %>%
        na.omit()
      
      d_sample2 <- d_sample[, c("sample_id", "band", "year", "min_age", "cap_doy", "age_group", "treatment",
                                "group", "mass", "bhead", "fwing", "b_cort", "s_cort", "d_cort",
                                "group_id", "comp_grp", "trt_label", "cross_grp")]
      
      pct_ad2 <- plyr::join(pct_ad2, d_sample2, "sample_id")
      pct_ad2$cpg_site <- paste(pct_ad2$chr, pct_ad2$start, pct_ad2$end, sep = "_")
      
      pct_be2 <- plyr::join(pct_be2, d_sample2, "sample_id")
      pct_be2$cpg_site <- paste(pct_be2$chr, pct_be2$start, pct_be2$end, sep = "_")
      
      
      bve_set <- unique(pct_be2$cpg_site)
      avd_set <- unique(pct_ad2$cpg_site)
      
      be <- filter(pct_be2, cpg_site %in% avd_set)
      ad <- filter(pct_ad2, cpg_site %in% bve_set)
      
      comb_pct <- rbind(be, ad)
      
      comb_pct2 <- comb_pct %>%
        pivot_wider(id_cols = c(band, treatment, cpg_site), names_from = group, values_from = value) %>%
        as.data.frame() %>%
        na.omit()
      
      colnames(comb_pct2) <- c("band", "treatment", "cpg_site", "b_post", "a_pre")
      
      site_list <- data.frame(cpg_site = unique(comb_pct2$cpg_site))
      site_list$rsq <- NA
      site_list$mod_sq <- NA
      site_list$err_sq <- NA
      
      
      for(i in 1:nrow(site_list)){
        sub <- subset(comb_pct2, comb_pct2$cpg_site == site_list$cpg_site[i])
        m <- lm(sub$b_post ~ sub$a_pre)
        site_list$rsq[i] <- summary(m)$r.squared
        aa <- anova(m)
        site_list$mod_sq[i] <- aa$`Mean Sq`[1]
        site_list$err_sq[i] <- aa$`Mean Sq`[2]
        print(i)
      }
      
      site_list$repeatability <- ((site_list$mod_sq - site_list$err_sq) / 2) / 
        (((site_list$mod_sq - site_list$err_sq) / 2) + site_list$err_sq)
      
      ggplot(sub, mapping = aes(x = a_pre, y = b_post)) +
        geom_point(mapping = aes(color = treatment)) +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
        xlab("Percent methylation pre-treatment") +
        ylab("Percent methylation post-treatment") +
        scale_color_manual(values = c("#1B9E77", "#D95F02")) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")
        
      
      
      cc <- comb_pct2 %>%
        group_by(cpg_site, treatment) %>%
        summarise(apre = mean(a_pre), bpost = mean(bpost))
      
      ggplot(cc, mapping = aes(x = apre, y = bpost)) +
        geom_hex() +
        #scale_color_viridis() +
        scale_fill_viridis() +
        #geom_smooth(method = "lm", inherit.aes = FALSE, mapping = aes(x = apre, y = bpost, color = treatment), se = FALSE)
        facet_wrap(~ treatment) +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
          xlab("Percent methylation pre-treatment") +
          ylab("Percent methylation post-treatment")