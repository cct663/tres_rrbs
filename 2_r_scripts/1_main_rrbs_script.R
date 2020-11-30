# Script for tree swallow RRBS data from Ithaca NY
# Written by Conor Taff: cct663@gmail.com ~ cct63@cornell.edu
# Last updated 11/2020

# Largely based on vignette here: http://127.0.0.1:18146/library/methylKit/doc/methylKit.html

# Load libraries ----
    pacman::p_load("methylKit", "tidyverse", "here")

# Load data ----
    cov_list <- list.files(here("0_processed_data/bismark_cov_output"))
    
# Pull sample id ----
    # merge to file that has filename and sample id
    
    
# Read into methylkit ----
    # Make the full file list
          file_list <- as.list(here("0_processed_data/bismark_cov_output", cov_list))
    
    # Read into methylrawlist
          meth_list <- methRead(file_list,
                                sample.id = as.list(cov_list),
                                assembly = "tres_20",
                                treatment = c(rep(c(1, 0), 60)),
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
          meth <- methylKit::unite(f_meth_list, destrand = FALSE, min.per.group = 10L)
          
    # Plot correlation in methylation between samples
          getCorrelation(meth, plot = TRUE)
          
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
          
          
          
          
          
          
          
          
    