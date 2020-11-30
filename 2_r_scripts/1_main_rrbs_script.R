# Script for tree swallow RRBS data from Ithaca NY
# Written by Conor Taff: cct663@gmail.com ~ cct63@cornell.edu
# Last updated 11/2020

# Largely based on vignette here: http://127.0.0.1:18146/library/methylKit/doc/methylKit.html

# Load libraries ----
    pacman::p_load("methylKit", "tidyverse", "here", "gridExtra", "ggpubr")

# Load data ----
    cov_list <- list.files(here("0_processed_data/bismark_cov_output"))
    d_bis <- read.delim(here("0_processed_data/bismark_sample_summary.txt"))
    d_sample <- read.delim(here("1_raw_data/rrbs_sample_metadata.txt"))
    
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
    
    
# Read into methylkit ----
    # Make the full file list
            file_list <- as.list(here("0_processed_data/bismark_cov_output", cov_list))
        
            test_list <- as.list(c(file_list[[1]], file_list[[25]], file_list[[49]],
                           file_list[[73]], file_list[[97]]))
            test_cov_list <- c(cov_list[1], cov_list[25], cov_list[49], cov_list[73], cov_list[97])
    
    # Read into methylrawlist
          meth_list <- methRead(file_list,
                                sample.id = as.list(cov_list),
                                assembly = "tres_20",
                                treatment = rep(c(1, 0), 60),
                                dbtype = "tabix",
                                context = "CpG",
                                pipeline = "bismarkCoverage",
                                header = FALSE,
                                dbdir = here("4_other_output"),
                                mincov = 10)
          
    # Descriptive stats of samples
        # Gives number of bases that fall into each bin of percent methylation
          # runs by sample. Could do as loop to make plots of all
            getMethylationStats(meth_list[[5]], plot = TRUE)
            
        # Gives coverage per sample. Same as above.
            getCoverageStats(meth_list[[5]], plot = TRUE)
            
    # Filters the list to exclude bases < x and > x coverage
        # greater than can be a problem because could be pcr bias
            f_meth_list <- filterByCoverage(meth_list, lo.count = 10, lo.perc = NULL,
                                            hi.count = NULL, hi.perc = 99.9)
            
    # Combine the whole methyllist together into one object
        # Set min per group to be minimum number of samples per group with coverage
          meth <- methylKit::unite(f_meth_list, destrand = FALSE, min.per.group = 1L)
          
    # Plot correlation in methylation between samples
          getCorrelation(meth, plot = TRUE)
          xx <- getCorrelation(meth)
          
    # Cluster samples
        #plot clustered dendrogram
          clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
          
        # PCA scree plot
          PCASamples(meth, screeplot = TRUE)
          PCASamples(meth)

    # Section on batch effects to look at (and potentially remove) PCs associated with covariates
          #sampleAnnotation
          #assocComp
          #removeComp
          #percMethylation
          #reconstruct
          
    # Tiling window analysis
        # tile the genome to 1k bp windows and step across, object can be fed into unite and 
            # calculateDiffMeth functions consecutively to get differential regions
            # may want to use a lower coverage limit for this
            
              # read in again with lower minimum
          
            # tile
              tiles = tileMethylCounts(meth, win.size = 1e4, step.size = 1e4, cov.bases = 10)
          
    # Finding differentially methylated bases or regions
        
        # Can do same thing here but feed in tiling output      
          myDiff <- calculateDiffMeth(meth)
            # can also turn on overdispersion = "MN" here
        
        # save hyper methylated regions   (change to hypo for opposite)
          diff25p_hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
         
        # Plot them, or turn off plot for text output 
          diffMethPerChr(myDiff, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25)
          
    # Next section details how to put in covariates to the model to check for differential
          
    # Next section uses library(genomation) to see where DMRs are compared to annotated genome
          #e.g., what percent are in promoters, etc
          
    # Most of remainder of code is little tweaks and file type conversions, etc.
          
          
          
          
          
          
          
          
    