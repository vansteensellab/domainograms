#!/usr/bin/Rscript

# Plot Domainograms (version 0.0.9)
## Vinicius Franceschini-Santos, Lise Dauban 230523

# - Functionality:
# This script generates the domainograms automatically from the GATC files. It 
# is standalone, so it contains all the LAD functions itselft. It plots all the
# standard domainograms (big, middle, and small) and also the scatterplots of the
# original Rmd file. You also can tell which allele the input files are from, so
# the domainogram will change -- considering that CAST allele contains a recom-
# bination. You can also set the type of recombination, as well as the start and
# end of it to get a proper domainogram.


# - Input
# (required)
# * GATC files (must be gziped)
# * annotation file
# * output directory
# * allele name (CAST or 129)
# * recombination type (will be considered only if --allele CAST)

# - Output (within output_dir)
# * scatterplots - raw and after loess correction
# * domainograms - big, middle, and small

# - Version
# 0.0.0 - First version
# 0.0.1 - Translate all plots to ggplot grammar
# 0.0.2 - Make damID window size as the first domainogram window size
# 0.0.3 - Implement ECD in the emp.domainogram function
# 0.0.4 - Fixed issue with domainogram plotting and missing data
# 0.0.5 - Fix gradient colors of domainogram
# 0.0.6 - Rework of the argument's name
# 0.0.7 - Handle human data
# 0.0.8 - Allow combined recombinations (deletion and inversion in same sample)
# 0.0.9 - Treat the control as inverted::experiment (for combined recombinations)
# 0.0.10 -Scale the logratios, converting to Z-scores

# Note: I had to load the LAD functions inside the main function. Maybe because
# of some trouble with similar object names or something like this. I decided not
# to dig into this, as load 

script_version <- "0.0.10"




set.seed(1818)
rlang::global_handle()
options(warn=-1, messages=-1)
suppressPackageStartupMessages({
        library("optparse")
        library(GenomicRanges)
        library(data.table) #for fast reading of data files
        library(rtracklayer)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(gggenes)
        library(patchwork)
        library(ggnewscale)
        library(stringr)
})
if (!require("ggbraid", quietly = TRUE)){
        cat("--------------------------- Installing ggbraid\n\n")
        remotes::install_github("nsgrantham/ggbraid")
} else {
        suppressPackageStartupMessages(library(ggbraid))
}


##############################################################################
## Main function #############################################################
##############################################################################

#parameters:

main = function(smooth_window_size, my_genes, del_start, del_end, inv_start, inv_end,
                annotation, exp_dam_files, exp_antibody_files, ctrl_dam_files, 
                ctrl_antibody_files, recombination_chr, output_dir, plot_title,
                script_version, del, inv, organism,
                recombination_centered, expression_tbl, M_L, M_R){
        
        # Creating output dir
        dir.create(output_dir, showWarnings = F)
        
        # logger <- function(s) {
        #         # Create logging messages
        #         cat(paste0(s, "\n"), file = paste0(output_dir, "/log_file.txt"), 
        #             append = T)
        #         cat(paste0(s, "\n"), file =  stderr())
        # }

        # COLORED LOGGER -- from https://github.com/r-lib/testthat/blob/717b02164def5c1f027d3a20b889dae35428b6d7/R/colour-text.r
        logger <- function(text, fg = "no", bg = NULL) {
  term <- Sys.getenv()["TERM"]
  colour_terms <- c("xterm-color","xterm-256color", "screen", "screen-256color")

  if(rcmd_running() || !any(term %in% colour_terms, na.rm = TRUE)) {
    cat(paste0(text), "\n", file =  stderr())
    cat(paste0(text, "\n"), file = paste0(output_dir, "/log_file.txt"), append=T)
  }

  col_escape <- function(col) {
    paste0("\033[", col, "m")
  }

  col <- .fg_colours[tolower(fg)]
  if (!is.null(bg)) {
    col <- paste0(col, .bg_colours[tolower(bg)], sep = ";")
  }

  init <- col_escape(col)
  reset <- col_escape("0")
  cat(paste0(init, text, reset, "\n"), file =  stderr())
  cat(paste0(text, "\n"), file = paste0(output_dir, "/log_file.txt"), append=T)

}

.fg_colours <- c(
  "black" = "0;30",
  "blue" = "0;34",
  "green" = "0;32",
  "cyan" = "0;36",
  "red" = "0;31",
  "purple" = "0;35",
  "brown" = "0;33",
  "light gray" = "0;37",
  "dark gray" = "1;30",
  "light blue" = "1;34",
  "light green" = "1;32",
  "light cyan" = "1;36",
  "light red" = "1;31",
  "light purple" = "1;35",
  "yellow" = "1;33",
  "white" = "1;37",
  "no" = "0"
)

.bg_colours <- c(
  "black" = "40",
  "red" = "41",
  "green" = "42",
  "brown" = "43",
  "blue" = "44",
  "purple" = "45",
  "cyan" = "46",
  "light gray" = "47"
)

rcmd_running <- function() {
  nchar(Sys.getenv('R_TESTS')) != 0
}
        
        logger("-----------------------------------------------------------", fg="green")
        logger(paste0("Plot domainograms wrapper script\nversion ", script_version), fg="green")
        logger("-----------------------------------------------------------", fg="green")
        
        logger("Command line:")
        logger(paste0("Rscript plot_domainograms.R",
                      "\n\t--exp_dam_files ", exp_dam_files,
                      "\n\t--exp_antibody_files ", exp_antibody_files,
                      "\n\t--ctrl_dam_files ", ctrl_dam_files,
                      "\n\t--ctrl_antibody_files ", ctrl_antibody_files,
                      "\n\t--annotation ", annotation,
                      "\n\t--output_dir ", output_dir,
                      "\n\t--left_margin ", M_L,
                      "\n\t--right_margin ", M_R,
                      "\n\t--plot_title \"", plot_title,
                      "\"\n\t--smooth_window_size ", smooth_window_size,
                      "\n\t--my_genes ", my_genes,
                      "\n\t--recombination_chr ", recombination_chr,
                      "\n\t--inv_start ", inv_start,
                      "\n\t--inv_end ", inv_end,
                      "\n\t--inv ", inv,
                      "\n\t--del_start ", del_start,
                      "\n\t--del_end ", del_end,
                      "\n\t--del ", del,
                      "\n\t--recombination_centered ", recombination_centered,
                      "\n\t--organism ", organism,
                      "\n\t--expression_tbl ", expression_tbl))
        logger("-----------------------------------------------------------")
        
        
        ########################################################################
        ## LAD functions -- sorry for building them inside the main :( 
        ##                           read the prologue of this script for details
        ########################################################################
        
        #superfast running sum calculator (from Elzo de Wit's domainogram function):
        #x is numeric vector
        #n is running window size
        runsum2<-function(x,n=21){
                cumsum(x)->sum.v
                sum.v<-c(0,sum.v)
                diff(sum.v,n)
        } #note: discards head and tail of size floor(n/2)!
        
        
        #1. Function to read all files and simply add read counts of replicates 
        #   for DamLam and Dam, respectively. Return as GRanges object
        read.combine.DamID.replicates<-function(DamLam, Dam){
                #DamLam (char): vector of file paths to DamLam replicates
                #Dam (char): vector of file paths to matching Dam-only replicates
        
                if(length(DamLam) != length(Dam)){
                        warning("warning: DamLam and Dam do not have the same\
                                 number of replicates")
                }
                #combine reads from DamLam datasets:
                for (i in 1:length(DamLam)){
                        data.df<-fread(cmd=paste("gunzip -c", DamLam[i]))
                        names(data.df)<-c("seqnames","start","end","DamLam")
                        logger(paste("DamLam replicate", i, "total read count:", 
                                      sum(data.df[,4])))
                        if(i==1){DL<-data.df$DamLam}
                        else{DL<-DL+data.df$DamLam}
                }
                for (i in 1:length(Dam)){
                        data.df<-fread(cmd=paste("gunzip -c", Dam[i]))
                        names(data.df)<-c("seqnames","start","end","Dam")
                        logger(paste("Dam replicate", i, "total read count:",
                                     sum(data.df[,4])))
                        if(i==1){D<-data.df$Dam}
                        else{D<-D+data.df$Dam}
                }
                data.df$Dam<-D; data.df$DamLam<-DL
                logger("converting to GRanges object...")
                GR<-makeGRangesFromDataFrame(data.df,keep.extra.columns = T)
                logger(paste("DamID data loaded and combined:", length(GR), 
                             "GATC fragments"))
                return(GR)
        }
        
        #2. Function to convert output of 1. into smoothed DamLam/Dam ratios. 
        #   Return as GRanges object
        smooth.DamID<-function(DamID.GR, title=NULL, ws=201, pseud=10)
                #returns GRanges object with smoothed log2(DamLam/Dam) as scores
                #note: it first does smoothing on Dam and DamLam each, then adds 
                #      pseudocount to each, THEN calculates log2-ratio
                #DamID.GR: GRanges object with two metadata columns: Dam and DamLam. 
                #          Generated by function read.combine.DamID.replicates()
                #ws (integer): window size for smoothing
                #pseud (integer): pseudocount to be added to each runsum2 
                #                 window prior to calculating the DamLam/Dam ratio 
                #                 in that window
        {
                stopifnot(ws%%2==1) #accept only odd window sizes
                #apply running window smoothing to Dam and DamLam each:
                Td<-runsum2(mcols(DamID.GR)$Dam, ws)
                Tdl<-runsum2(mcols(DamID.GR)$DamLam, ws) 
                #add padding at beginning and add to compensate for lost data
                pad<-floor(ws/2)
                Td<-c(rep(NA,pad), Td, rep(NA,pad))
                Tdl<-c(rep(NA,pad), Tdl, rep(NA,pad))
                #set elements to NA that have no reads in either Dam or DamLam:
                nd<-which(Td==0 & Tdl==0)
                Td[nd]<-NA; Tdl[nd]<-NA
                #add pseudocount to both channels: 
                Td<-Td+pseud
                Tdl<-Tdl+pseud
                #calculate logratios, normalize to mean:
                Tr<-log2(Tdl/Td)
                Tr<-Tr-mean(Tr, na.rm=TRUE)
                mcols(DamID.GR)$score<-Tr
                #mcols(DamID.GR)<-mcols(DamID.GR)[,3]
                if(is.null(title)){
                       #mcols(DamID.GR)<-mcols(DamID.GR)[,3] 
                       mcols(DamID.GR)<-as.numeric(scale(mcols(DamID.GR)[,3]))
                } else {
                        logger(paste("1 SD for the", title, "is", sd(mcols(DamID.GR)[,3], na.rm=T), "(log2ratio units)"))
                        mcols(DamID.GR)<-as.numeric(scale(mcols(DamID.GR)[,3]))                                                  #################################### add scaling
                }
                names(mcols(DamID.GR))<-"logratio"
                return(DamID.GR) #returns Granges object with smoothed log2(DamLam/Lam) 
                #                 as metadata
        }
        
        
        #3A. Domainogram calculation based on output .1
        #note: This may not be statically correct, because the smoothed windows overlap
        #in the original de Wit et al PLoS Genet 2008 paper this was not the case.
        #Use at own risk! 
        #See emp.domainogram() function below; which is an alternative to simply 
        # calculate genome-wide rankings of domain scores. 
        #DOMAINOGRAM function:
        domainogram.bvs <- function(damid.exp, damid.ctl, chrom, 
                                    window.sizes=c(101), pseud=10,
                                    up=FALSE, logratios=FALSE)
                #input: 
                #damid.exp and damid.ctl: GRanges objects with DamID data as 
                #                         generated by read.combine.damid.replicates()
                #chrom: chromosome to be analyzed
                #window.sizes: vector of odd integers, window.sizes to be 
                #              used for domainogram
                #pseud: pseudocount used in calculation of logratios
                #up: should be tested for signals in Exp going up compared 
                #    to Ctl? Set to FALSE if testing for signals going down.
                #logratios: should calculations be done on DamID logratios
                #           (TRUE) or converted to linear ratios (FALSE)?
                #output: GRanges object with domainogram data.frame as metadata
                #        (accessible via mcols() function)
        {
                stopifnot(sum(window.sizes%%2)==length(window.sizes))
                for(w.size in window.sizes){
                        cat("calculating window size", w.size, "\r")
                        #first calculate (log)ratios using current window size
                        Ctl<-smooth.DamID(damid.ctl, ws=w.size, pseud=pseud, title="control")
                        Exp<-smooth.DamID(damid.exp, ws=w.size, pseud=pseud, title="experiment")
                        #extract data vector:
                        Ctl<-mcols(Ctl)$logratio
                        Exp<-mcols(Exp)$logratio
                        #convert to linspace if logratios==FALSE:
                        if(logratios==FALSE) {Ctl<-2^Ctl; Exp<-2^Exp}
                        #calculate difference and convert to quantiles:
                        diff.damid <- Exp-Ctl
                        if(up==FALSE) {diff.damid<- -diff.damid}
                        q.vec <- (rank(diff.damid)-0.5)/length(diff.damid)
                        #calculate p.values:
                        run.q <- runsum2(log(q.vec[as.character(seqnames(damid.exp))==chrom]),
                                         n=w.size )
                        chi.val <- -2*run.q
                        comb.p <- pchisq(chi.val, df=w.size*2)
                        comb.p <- -log10(comb.p)
                        #make sure that comb.p is the same length again as 
                        #input score vector:
                        pad<-floor(w.size/2)
                        comb.p<-c(rep(NA,pad), comb.p, rep(NA,pad))
                        #put results into a matrix:
                        if(w.size==window.sizes[1]){
                                p.matrix<-comb.p
                        } else {
                                p.matrix<-cbind(p.matrix, comb.p)
                        }
                } #end loop through window.sizes
                
                #return results as GRanges object:
                GR<-damid.ctl[seqnames(damid.ctl)==chrom] #create GRanges object for chrom only
                mcols(GR)<-p.matrix #add domainogram data as metadata
                names(mcols(GR))<-as.character(window.sizes)
                return(GR) 
        }
        
        #3B. Empirical domainogram calculation based on output .1
        #this differs from the regular domainogram in that the 
        # null distribution is entirely empirical, 
        #based on the entire genome.
        #the returned values are simply the ranks of the differences
        # in DamID values, for windows of the same size & genome-wide
        emp.domainogram <- function(damid.exp, damid.ctl, chrom, window.sizes=c(101), 
                                    pseud=10, up=FALSE, logratios=FALSE, 
                                    loess.corr=NULL, baseline=0)
                #damid.exp and damid.ctl: GRanges objects with DamID data as generated
                # by read.combine.damid.replicates()
                #chrom: chromosome to be analyzed
                #window.sizes: vector of odd integers, window.sizes to be
                # used for domainogram
                #pseud: pseudocount used in calculation of logratios
                #up: should be tested for signals in Exp going up compared to Ctl?
                # Set to FALSE if testing for signals going down.
                #logratios: should calculations be done on DamID logratios (TRUE)
                # or converted to linear ratios (FALSE)?
                #loess.corr: a loess object to straighten the overall exp vs
                # control relationship
                #baseline: if logratios==TRUE, the quantile below which data
                # are not included in the rank scores.
                #-- this is because changes in DamID data for very low values
                # (in log space) are essentially baseline fluctuations that
                # are not relevant
                #output: GRanges object with domainogram data.frame as
                # metadata (accessible via mcols() function)
        {
                stopifnot(sum(window.sizes%%2)==length(window.sizes))
                for(w.size in window.sizes){
                        logger(paste0("Working with window ", w.size))
                        #first calculate (log)ratios using current window size
                        Ctl<-smooth.DamID(damid.ctl, ws=w.size, pseud=pseud)
                        Exp<-smooth.DamID(damid.exp, ws=w.size, pseud=pseud)
                        #extract data vector:
                        Ctl<-mcols(Ctl)$logratio
                        Exp<-mcols(Exp)$logratio
                        #apply loess correction:
                        if(!is.null(loess.corr)) {Exp <- Exp - predict(loess.corr, Ctl)}
                        #convert to linspace if logratios==FALSE:
                        if(logratios==FALSE) {Ctl<-2^Ctl; Exp<-2^Exp}
                        #calculate difference :
                        diff.damid <- Exp-Ctl
                        if(up==FALSE) {diff.damid<- -diff.damid}
                        #mask baseline effects if working in log space:
                        #baseline is defined as all values below 0.3 quantile
                        if(logratios==TRUE) {
                                diff.damid<-ifelse(
                                        Exp< quantile(Exp, baseline, na.rm=TRUE) & 
                                        Ctl< quantile(Ctl, baseline, na.rm=TRUE), 
                                        NA, 
                                        diff.damid)}
                        # #convert to quantiles
                        # vf20221201: changed the old calculation to actual ECDF.
                        #             It's a more robust way of ranking the data
                        q.vec = ecdf(diff.damid)(diff.damid)
                        # q.vec <- (rank(diff.damid, na.last="keep")
                        #           /
                        #           length(diff.damid[!is.na(diff.damid)]))
                        
                        # #calculate q.values:
                        #put results into a matrix:
                        if(w.size==window.sizes[1]){
                                q.matrix<-q.vec[as.character(seqnames(damid.exp))==chrom]
                        } else {
                                q.matrix<-cbind(
                                        q.matrix,
                                        q.vec[as.character(seqnames(damid.exp))==chrom])
                        }
                        logger("-- Done")
                } #end loop through window.sizes
                #return results as GRanges object:
                GR<-damid.ctl[seqnames(damid.ctl)==chrom] #create GRanges object for chrom only
                mcols(GR)<-q.matrix #add domainogram data as metadata
                names(mcols(GR))<-as.character(window.sizes)
                return(GR) 
        }
        
        
        
        #4. Function to plot domainorams with ggplot
        #expression_tbl = "~/mydata/GitHub/plot_domainograms/4DNFIPFKK5LM_RNAseq_fixedIDs.tsv"

        plot_domainograms_with_ggplot = function(dg_exp, del, inv, annot, my_genes, 
                                                 recombination_centered, mrg.left,
                                                 mrg.right, smooth_window_size,
                                                 sDamID.Exp, sDamID.Control, logratio=T,
                                                 y_limit=NA, expression_tbl, organism,
                                                 plot_title, Del, Inv, inv_end, inv_start, del_end, del_start){
                # recombination_centered: a TRUE/FALSE object indicating if the plot 
                #                         should be centered in the recombination or not.
                #                         If FALSE, genes in 'my_genes' will be the center
                ## Determine gene and plotting coordenates
                rows_of_your_genes = which(annot$gene %in% my_genes)
                chr_of_your_genes = unique(annot$chrom[rows_of_your_genes])
                # Stop if gene appears to be on multiple chromosomes
                stopifnot(length(chr_of_your_genes)==1) 
                
                # Changing the x-range of the plot according to 'recombination_centered'
                #if(isTRUE(recombination_centered)){
                #        x_range<-c(recombination_start - mrg.left,
                #                   recombination_end + mrg.right)
                #} else {
                #        x_range<-c(min(annot$start[rows_of_your_genes]) - mrg.left,
                #                   max(annot$stop[rows_of_your_genes]) + mrg.right)
                #}
                x_range<-c(min(annot$start[rows_of_your_genes]) - mrg.left,
                                   max(annot$stop[rows_of_your_genes]) + mrg.right)
                # Adjusting the x-range
                if(x_range[1]<1) {x_range[1]<-1}
                if(x_range[2]>max(end(dg_exp))) {x_range[2]<-max(end(dg_exp))}
                
                genes_in_your_chr = annot[annot$chrom==chr_of_your_genes,] 
                
                # Transforming the smooth window size (swindorsize) to kb. (previously,
                #   it was expressed in number of GATC fragments)
                swinsizekb = 
                        # fragments within the x_range
                        dg_exp[start(dg_exp)>=x_range[1] & end(dg_exp)<=x_range[2]] %>% 
                        # get their size
                        width() %>% 
                        # get mean of their size, and transform to kb
                        mean() * smooth_window_size / 1000
                # round the result
                swinsizekb = round(swinsizekb,0)
                
                ########################################################################
                #### PLOT PADAMID PAIR
                
                # get genome-wide quantiles
                Q005 = quantile(
                        (mcols(sDamID.Exp)$logratio 
                         + mcols(sDamID.Control)$logratio) / 2,
                        0.05,
                        na.rm=TRUE)
                
                Q095 = quantile(
                        (mcols(sDamID.Exp)$logratio
                         + mcols(sDamID.Control)$logratio) / 2,
                        0.95,
                        na.rm=TRUE)
                
                # extract relevant part of the data:
                neg = sDamID.Exp                                        ##############
                neg$logratio = -3                                       ##############
                if(del == TRUE){
                        # construct data frame for deletions
                        # Explanations on this files can be found after this `if`
                        Exp<-sDamID.Exp[(seqnames(sDamID.Exp)==chr_of_your_genes
                                         & start(sDamID.Exp)>=x_range[1] 
                                         & end(sDamID.Exp)<=x_range[2])]
                        
                        Exp_before_rec <-sDamID.Exp[(seqnames(sDamID.Exp)==chr_of_your_genes
                                                     & start(sDamID.Exp)>=x_range[1] 
                                                     & end(sDamID.Exp)<del_start)]
                        
                        Exp_after_rec <-sDamID.Exp[(seqnames(sDamID.Exp)==chr_of_your_genes
                                                    & start(sDamID.Exp)>del_end 
                                                    & end(sDamID.Exp)<=x_range[2]
                                                    & !is.na(mcols(sDamID.Exp)$logratio))]
                        
                        Rec<-neg[(seqnames(neg)==chr_of_your_genes 
                                  & start(neg)>=del_start 
                                  & end(neg)<=del_end)]
                        
                        
                        
                        Ctl<-sDamID.Control[(seqnames(sDamID.Control)==chr_of_your_genes
                                             & start(sDamID.Control)>=x_range[1]
                                             & end(sDamID.Control)<=x_range[2])]
                        
                        
                        damid_df = rbind(
                                data.frame(chr_position = start(Exp),
                                           damid_score = mcols(Exp)$logratio,
                                           library = "Experiment"),
                                
                                data.frame(chr_position = start(Ctl),
                                           damid_score = mcols(Ctl)$logratio,
                                           library = "Control")
                        )
                        
                        damid_df_b = rbind(
                                data.frame(chr_position = start(Exp_before_rec),
                                           damid_score = Exp_before_rec$logratio,
                                           library = "Experiment"),
                                
                                data.frame(chr_position = start(Ctl),
                                           damid_score = mcols(Ctl)$logratio,
                                           library = "Control")
                        )
                        
                        damid_df_a = rbind(
                                data.frame(chr_position = start(Exp_after_rec),
                                           damid_score = Exp_after_rec$logratio,
                                           library = "Experiment"),
                                
                                data.frame(chr_position = start(Ctl),
                                           damid_score = mcols(Ctl)$logratio,
                                           library = "Control")
                        )
                        
                        dealing_with_deletions = function(x){
                                if(length(x)>1){
                                        unlist(x)[2]
                                }else{
                                        unlist(x)[1]
                                }
                                
                        }
                        
                        damid_df_wide =  pivot_wider(damid_df, 
                                                     names_from = library, 
                                                     values_from = damid_score,
                                                     values_fn = dealing_with_deletions)
                        
                        damid_df_wide_b =
                                damid_df_wide %>%
                                filter(chr_position <= del_start)
                        
                        damid_df_wide_a =
                                damid_df_wide %>%
                                filter(chr_position >= del_end)
                        
                        damid_df_fill_b = 
                                damid_df_wide_b %>%
                                mutate(ymax = pmin(Control, Experiment, na.rm = T))
                        
                        damid_df_fill_a = 
                                damid_df_wide_a %>%
                                mutate(ymax = pmin(Control, Experiment, na.rm = T))
                        
                        
                        
                        
                } else {
                        # construct data frame for inversions and non-recombinant data
                        # Explanations on this files can be found after this `if`
                        Exp<-sDamID.Exp[(seqnames(sDamID.Exp)==chr_of_your_genes
                                         & start(sDamID.Exp)>=x_range[1] 
                                         & end(sDamID.Exp)<=x_range[2])]
                        
                        Ctl<-sDamID.Control[(seqnames(sDamID.Control)==chr_of_your_genes
                                             & start(sDamID.Control)>=x_range[1]
                                             & end(sDamID.Control)<=x_range[2])]
                        
                        Neg<-neg[(seqnames(neg)==chr_of_your_genes 
                                  & start(neg)>=x_range[1] 
                                  & end(neg)<=x_range[2])]
                        
                        damid_df = rbind(
                                data.frame(chr_position = start(Exp),
                                           damid_score = mcols(Exp)$logratio,
                                           library = "Experiment"),
                                
                                data.frame(chr_position = start(Ctl),
                                           damid_score = mcols(Ctl)$logratio,
                                           library = "Control")
                        )
                        
                        damid_df_wide =  pivot_wider(damid_df, 
                                                     names_from = library, 
                                                     values_from = damid_score)
                        damid_df_fill = 
                                damid_df_wide %>%
                                mutate(ymax = pmin(Control, Experiment, na.rm = T))        
                }
                
                
                
                # Adjust if not logratio:
                if(logratio==FALSE){
                        mcols(Exp)$logratio<-2^mcols(Exp)$logratio
                        mcols(Ctl)$logratio<-2^mcols(Ctl)$logratio
                        Q005<-2^Q005; Q095<-2^Q095
                }
                
                # Adjust if y_limit is not set
                if(is.na(y_limit)){
                        y_limit<-range(Q005, Q095, mcols(Exp)$logratio,
                                       mcols(Ctl)$logratio, na.rm=TRUE)
                        if(logratio==FALSE){y_limit[1]<-0}
                }
                
                
                ########################################################################
                ## (1.1) Plotting DamID
                ########################################################################
                #### (a) Contrustucting dataframes
                ########################################################################
                
                ## overview of damid_df:
                #
                #       chr_position  damid_score   library
                # 1        108270273  -1.516950  Experiment
                # 2        108270320  -1.516950  Experiment
                # 3        108270419  -1.516950  Experiment
                # 4        108270707  -1.516950  Experiment
                # 5        108270879  -1.493826  Experiment
                # 6        108271737  -1.492966  Experiment
                # .        ...        ...         ...
                # 28437    113923706  -0.8189863 Control
                # 28438    113924544  -0.7992707 Control
                # 28439    113925351  -0.7864777 Control
                # 28440    113925812  -0.8113435 Control
                # 28441    113925974  -0.8113435 Control
                # 28442    113926789  -0.8263779 Control
                
                
                ## overview of damid_df_wide:
                #  
                #
                #   chr_position Experiment Control
                #        <int>      <dbl>     <dbl>
                # 1    108270273      -1.52   -1.68
                # 2    108270320      -1.52   -1.68
                # 3    108270419      -1.52   -1.68
                # 4    108270707      -1.52   -1.68
                # 5    108270879      -1.49   -1.66
                # 6    108271737      -1.49   -1.66
                # .      ...             .       .
                # 1    113923706     -0.482  -0.819
                # 2    113924544     -0.462  -0.799
                # 3    113925351     -0.450  -0.786
                # 4    113925812     -0.461  -0.811
                # 5    113925974     -0.461  -0.811
                # 6    113926789     -0.473  -0.826
                
                
                
                ## overview of damid_df_fill:
                #
                #
                # chr_position Experiment Control  ymax
                # <int>      <dbl>   <dbl> <dbl>
                # 1    108270273      -1.52   -1.68 -1.68
                # 2    108270320      -1.52   -1.68 -1.68
                # 3    108270419      -1.52   -1.68 -1.68
                # 4    108270707      -1.52   -1.68 -1.68
                # 5    108270879      -1.49   -1.66 -1.66
                # 6    108271737      -1.49   -1.66 -1.66
                # .       ...           .       .    .
                # 1    113923706     -0.482  -0.819 -0.819
                # 2    113924544     -0.462  -0.799 -0.799
                # 3    113925351     -0.450  -0.786 -0.786
                # 4    113925812     -0.461  -0.811 -0.811
                # 5    113925974     -0.461  -0.811 -0.811
                # 6    113926789     -0.473  -0.826 -0.826
                
                
                ## overview of segment:
                #
                # - its just a tricky to plot the window size inside the chart
                # x        y      xend     yend
                # 95% 108976958 1.387198 109097958 1.387198
                segment = data.frame(x = diff(x_range)/8+x_range[1], 
                                     y = Q095-0.1,
                                     xend = diff(x_range)/8+x_range[1]+(swinsizekb*1e3), 
                                     yend = Q095-0.1)
                
                #### (b) plotting
                ########################################################################
                
                if(del == TRUE){
                        damid_plot =
                                ggplot() +
                                # First braid: filling the area between Control and Experiment
                                geom_braid(aes(chr_position, 
                                               ymin = Control, 
                                               ymax = Experiment, 
                                               fill = Control < Experiment),
                                           alpha = 0.5,
                                           data = damid_df_wide_b,
                                           method = "line") +
                                geom_braid(aes(chr_position, 
                                               ymin = Control, 
                                               ymax = Experiment, 
                                               fill = Control < Experiment),
                                           alpha = 0.5,
                                           data = damid_df_wide_a,
                                           method = "line") +
                                # Second braid: filling the area bellow the lines
                                geom_braid(aes(chr_position, 
                                               ymin = y_limit[1],
                                               ymax = ymax), 
                                           fill = "#f2f2f2",
                                           data = damid_df_fill_b,
                                           method = "line") +
                                geom_braid(aes(chr_position, 
                                               ymin = y_limit[1],
                                               ymax = ymax), 
                                           fill = "#f2f2f2",
                                           data = damid_df_fill_a,
                                           method = "line") +
                                # Lines of DamID score 
                                geom_line(aes(chr_position, damid_score, color = library),
                                          alpha = 1,
                                          data = damid_df_b) +
                                geom_line(aes(chr_position, damid_score, color = library),
                                          alpha = 1,
                                          data = damid_df_a) +
                                geom_line(aes(x=x_range[1], y=0, alpha = "tmp"), size = 2) +
                                # Changing colors
                                scale_fill_manual(values = c("#3c5488", "#e64b35")) +
                                geom_segment(aes(x = x, 
                                                 y = y,
                                                 xend = xend, 
                                                 yend = yend),
                                             size = 2,
                                             alpha = 1,
                                             lineend = "butt",
                                             colour = "black",
                                             data = segment,
                                             show.legend = F) +
                                # Add dashed lines
                                geom_vline(xintercept = del_start,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = del_end,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = inv_start,
                                   linetype = "dashed",
                                   color = "#b09c85") +
                                geom_hline(yintercept = Q005,
                                           linetype = "dotted",
                                           color = "gray80") +
                                geom_hline(yintercept = Q095,
                                           linetype = "dotted",
                                           color = "gray80") +
                                ggthemes::theme_few() +
                                ylab("pA-DamID (z-score)") +
                                xlab(chr_of_your_genes)  +
                                scale_x_continuous(limits = x_range,
                                                   expand = c(0, 0)) +
                                
                                theme(axis.title.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.line.y.left = element_line(color = "black"),
                                      axis.line.y.right = element_line(color = "black")) +
                                
                                scale_y_continuous(limits = y_limit,
                                                   expand = c(0,0.1))+
                                
                                scale_color_manual(values = c("#3c5488", "#e64b35"),
                                                   name = "Library",
                                                   labels = c(data.summ[1], data.summ[2]),
                                                   guide = guide_legend(override.aes = list(width = 20,
                                                                                            alpha = 1)))+
                                scale_alpha_manual(values = 1,
                                                   name = NULL,
                                                   labels = "Sliding window size",
                                                   guide = guide_legend(override.aes = list(alpha = 1)))+
                                guides(fill = "none")
                } else {
                        damid_plot = 
                                ggplot() +
                                # First braid: filling the area between Control and Experiment
                                geom_braid(aes(chr_position, 
                                               ymin = Control, 
                                               ymax = Experiment, 
                                               fill = Control < Experiment),
                                           alpha = 0.5,
                                           data = damid_df_wide,
                                           method = "line") +
                                # Second braid: filling the area bellow the lines
                                geom_braid(aes(chr_position, 
                                               ymin = y_limit[1],
                                               ymax = ymax), 
                                           fill = "#f2f2f2",
                                           data = damid_df_fill,
                                           method = "line") +
                                # Lines of DamID score 
                                geom_line(aes(chr_position, damid_score, color = library),
                                          alpha = 1,
                                          data = damid_df) +
                                geom_line(aes(x=x_range[1], y=0, alpha = "tmp"), size = 2) +
                                # Changing colors
                                scale_fill_manual(values = c("#3c5488", "#e64b35")) +
                                geom_segment(aes(x = x, 
                                                 y = y,
                                                 xend = xend, 
                                                 yend = yend),
                                             size = 2,
                                             alpha = 1,
                                             lineend = "butt",
                                             colour = "black",
                                             data = segment,
                                             show.legend = F) +
                                # Add dashed lines
                                geom_vline(xintercept = del_start,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = del_end,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = inv_start,
                                   linetype = "dashed",
                                   color = "#b09c85") +
                                geom_hline(yintercept = Q005,
                                           linetype = "dotted",
                                           color = "gray80") +
                                geom_hline(yintercept = Q095,
                                           linetype = "dotted",
                                           color = "gray80") +
                                ggthemes::theme_few() +
                                ylab("pA-DamID (z-score)") +
                                xlab(chr_of_your_genes)  +
                                scale_x_continuous(limits = x_range,
                                                   expand = c(0, 0)) +
                                
                                theme(axis.title.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.line.y.left = element_line(color = "black"),
                                      axis.line.y.right = element_line(color = "black")) +
                                
                                scale_y_continuous(limits = y_limit,
                                                   expand = c(0,0.1))+
                                
                                scale_color_manual(values = c("#3c5488", "#e64b35"),
                                                   name = "Library",
                                                   labels = c(data.summ[1], data.summ[2]),
                                                   guide = guide_legend(override.aes = list(width = 20,
                                                                                            alpha = 1)))+
                                scale_alpha_manual(values = 1,
                                                   name = NULL,
                                                   labels = "Sliding window size",
                                                   guide = guide_legend(override.aes = list(alpha = 1)))+
                                guides(fill = "none")
                }
                
                ########################################################################
                ## (1.2) Plotting domainograms
                ########################################################################
                #### (a) constructing dataframes
                ########################################################################
                
                # subset of domainograms withing the region of interest
                dgsub = DG.Exp[(start(DG.Exp)>=x_range[1] 
                                & end(DG.Exp)<=x_range[2])]
                
                # transforming the output of 'calculate domainograms' in dataframe
                genome_wide_df = as.data.frame(mcols(DG.Exp))
                focused_df = as.data.frame(mcols(dgsub))
                
                # things for the y=axis
                gatc2len<-mean(width(dgsub)) #average gatc fragment length
                wsizes<-log10(as.integer(names(mcols(dgsub)))*gatc2len)
                
                # Coloring based on ECDF values.
                #  -- with this, the color of domainograms are a gradient ranging from 
                #     0.95 to 1 for the control data and from 0.05 to 0 for the 
                #     experimental data. This values mean which quantile each number 
                #     belongs to
                
                ## overview of ecdf_exp:
                #
                #
                # start        ECD Window
                # 644  108506204 0.03313090     67
                # 645  108506591 0.04997497     67
                # 1393 108793978 0.04781274     67
                # 1394 108795402 0.04781274     67
                # 1395 108796095 0.03505622     67
                # 1396 108798612 0.03505622     67
                # . ..   . . .     . .       .
                # 1     NA  NA   1451
                # 11    NA  NA   1669
                # 12    NA  NA   1919
                # 13    NA  NA   2207
                # 14    NA  NA   2537
                # 15    NA  NA   2917
                
                ## overview of ecdf_ctl
                #
                #
                # start       ECD Window
                # 416 108424466 0.9524242     67
                # 417 108424909 0.9508959     67
                # 422 108426177 0.9504722     67
                # 423 108426214 0.9504722     67
                # 424 108426351 0.9504722     67
                # 426 108426424 0.9507825     67
                # . ..   . . .     . .       .
                # 1     NA  NA   1451
                # 11    NA  NA   1669
                # 12    NA  NA   1919
                # 13    NA  NA   2207
                # 14    NA  NA   2537
                # 15    NA  NA   2917
                
                ecdf_exp = data.frame()
                ecdf_ctl = data.frame()
                for(this_col in 1:length(colnames(focused_df))){ 
                        #genome_wide_not2 = genome_wide_df[genome_wide_df!= 2,this_col]
                        
                        # estimating ECD function
                        # ecd_function = 
                        #         genome_wide_not2 %>%
                        #         na.omit() %>%
                        #         ecdf()
                        # 
                        # get start positions 
                        start_pos = start(dgsub)[mcols(dgsub)[,this_col]!=2]
                        # get domainogram scores 
                        scores = mcols(dgsub)[,this_col][mcols(dgsub)[,this_col]!=2]
                        
                        # convert window to int. first output is NA. second is the colname
                        window_name = 
                                colnames(focused_df)[this_col] %>%
                                strsplit("X") %>%
                                unlist() %>%
                                as.numeric()
                        
                        tmp = data.frame(start = start_pos, 
                                         # ecd of each score
                                         score = scores,
                                         #ECD = ecd_function(scores),
                                         # Window is one of 67, 77, ...
                                         Window = window_name[2])
                        
                        # logger(paste("DEBUG:: Colum", window_name[2]))
                        # logger("Debug:::ECDF dataframe - CTL")
                        # logger(tmp %>% arrange(ECD) %>% head() %>% knitr::kable(align = "c"))
                        # logger("Debug:::ECDF dataframe - EXP")
                        # logger(tmp %>% arrange(-ECD) %>% head() %>% knitr::kable(align = "c"))
                        # logger(paste("------------------------------"))
                        
                        # complete with NA if this window does not contain significant data
                        if(nrow(tmp[tmp$score>=0.95,]) == 0){
                                higher_score = tmp %>% arrange(-score) %>% head(n=1)
                                logger(paste0("Window ", 
                                             window_name[2], 
                                             "'s highest score is: ",
                                             higher_score$score%>% round(5) * 100,
                                             "% (below 95%)"), fg="yellow")
                                #logger(higher_score%>% knitr::kable(align = "c"))
                                tmp_df = data.frame(start = 0, 
                                                 #ECD = NA,
                                                 score = 0.5,
                                                 Window = window_name[2])
                                # control DF containing values in 5% quantile or less
                                ecdf_ctl = rbind(ecdf_ctl, tmp_df)
                                ecdf_exp = rbind(ecdf_exp, tmp[tmp$score<=0.05,])
                        #} else if(!any(tmp$ECD<=0.05)){
                         } else if(nrow(tmp[tmp$score<=0.05,]) == 0){
                                 lowest_score = tmp %>% arrange(score) %>% head(n=1)
                                 logger(paste0("Window ", 
                                               window_name[2], 
                                               "'s lowest score is: ",
                                               lowest_score$score %>% round(5)*100,
                                               "% (above 5%)"), fg="yellow")
                                tmp_df = data.frame(start = 0, 
                                                    # ecd of each score
                                                    #ECD = NA,
                                                    score=0.5,
                                                    # Window is one of 67, 77, ...
                                                    Window = window_name[2])
                                #logger(lowest_score %>% knitr::kable(align = "c"))
                                
                                # experimental DF containing values in the 5% or more quantile
                                ecdf_exp = rbind(ecdf_exp, tmp_df)
                                ecdf_ctl = rbind(ecdf_ctl, tmp[tmp$score>=0.95,])
                        } else {
                                logger(paste0("Window ", 
                                              window_name[2], 
                                              " contains scores within both quantiles"))
                                # experimental DF containing values in the 5% quantile or less
                                ecdf_exp = rbind(ecdf_exp, tmp[tmp$score<=0.05,])
                                # control DF containing values in 95% or more quantile
                                ecdf_ctl = rbind(ecdf_ctl, tmp[tmp$score>=0.95,]) 
                        }
                }
                
                # print the window 829 to logger ( to check the truncated data )
                # options(datatable.print.nrows=10)
                # logger("Debug:::ECDF dataframe - CTL")
                # logger(ecdf_ctl[ecdf_ctl$Window == 829,] %>% arrange(ECD) %>% head() %>% knitr::kable())
                # logger(ecdf_ctl[ecdf_ctl$Window == 829,] %>% arrange(ECD) %>% tail() %>% knitr::kable())
                # logger("Debug:::ECDF dataframe - EXP")
                # logger(ecdf_exp[ecdf_exp$Window == 829,] %>% arrange(ECD) %>% head() %>% knitr::kable())
                # logger(ecdf_exp[ecdf_exp$Window == 829,] %>% arrange(ECD) %>% tail() %>% knitr::kable())
                
                # create the recombinant dataframe
                if(inv == TRUE | del == TRUE){
                        recomb = data.frame()
                        for(this_col in 1:length(colnames(focused_df))){ 
                                # getting those values == 2 (Lise trick to deal 
                                # with recombinations)
                                recombinations = focused_df[,this_col] == 2
                                
                                # get start positions
                                start_pos = start(dgsub)[recombinations]
                                
                                # convert window to int. first output is NA. 
                                # second is the colname
                                window_name = 
                                        colnames(focused_df)[this_col] %>%
                                        strsplit("X") %>%
                                        unlist() %>%
                                        as.numeric()
                                
                                tmp = data.frame(start = start_pos,
                                                 # no need ECD 
                                                 Window = window_name[2])
                                recomb = rbind(recomb, tmp)
                        }
                }
                
                ## Tricky thing to get the right axis to show...
                tmp = round(wsizes,1)
                a=0
                for(i in 1:length(tmp)){if(tmp[i]%%0.5==0){a[i]=tmp[i]}else{a[i]=""}}
                
                b = a
                for(i in seq(1, length(b), 2)){if(is.na(b[i]) | is.na(b[i+1])) next
                        else if(b[i] == b[i+1]){b[i+1] = ""}
                }
                b[duplicated(b)] = ""
                
                #### (b) plotting
                ########################################################################
                # VF230804: assess both experiment and control windows to get the
                #           y-axis for the domainograms
                y_axis = c(ecdf_exp$Window, ecdf_ctl$Window) %>% factor() %>% levels()
                domainogram_plot = 
                        ggplot()+
                        geom_tile(aes(x = start, 
                                      y = factor(Window), 
                                      color = score), 
                                  data = ecdf_exp) +
                        scale_color_distiller(palette = "Reds", 
                                              limits = c(0,0.05)) +
                        new_scale_color()+
                        geom_tile(aes(x = start, 
                                      y = factor(Window), 
                                      color = score), 
                                  data = ecdf_ctl) +
                        scale_color_distiller(palette = "Blues", direction = 1,
                                              limits = c(0.95, 1)) +
                        ggthemes::theme_few()+
                        scale_x_continuous(limits = x_range,
                                           expand = c(0, 0)) +
                        scale_y_discrete(
                                # removes 'NA' from the axis
                                limits = y_axis,
                                # get this labels...
                                breaks = y_axis[b != ""],
                                # ...and rename it to this
                                labels = b[b != ""]) +
                        geom_vline(xintercept = del_start,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = del_end,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = inv_start,
                                   linetype = "dashed",
                                   color = "#b09c85") +
                        theme(axis.title.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              legend.position = "none",
                              axis.line.y.left = element_line(color = "black"),
                              plot.title = element_text(hjust = 0.5, size = 15)) +
                        ggtitle(plot_title) +
                        ylab(expression("log"[10]~"(scale [bp])"))
                
                if(inv == TRUE | del == TRUE){
                        domainogram_plot = 
                                domainogram_plot + 
                                geom_tile(aes(x = start, 
                                              y = factor(Window)), 
                                          data = recomb,
                                          color = "#f2f2f2")
                }
                ########################################################################
                ## (1.3) plotting genes
                ########################################################################
                #### (a) creating dataframes
                ########################################################################
                
                # Parsing annotation
                genes_to_plot = annot[(annot$start>=x_range[1] 
                                       & annot$stop<=x_range[2]
                                       & annot$chrom==chr_of_your_genes),]
                # Getting only the occurrence of genes (remove sub genes)
                genes_to_plot.first = genes_to_plot[match(unique(genes_to_plot$gene), 
                                                          genes_to_plot$gene),]
                # Get geneID for those in plot.first
                # this step is organism-specific, so it need to change if its human data
                if(organism == "mouse"){
                        suppressPackageStartupMessages(library(org.Mm.eg.db))
                        genes_to_plot.first$Ensembl = mapIds(org.Mm.eg.db, 
                                                             keys = genes_to_plot.first$gene, 
                                                             keytype = "SYMBOL", 
                                                             column = "ENSEMBL")
                        
                } else if (organism == "human") {
                        suppressPackageStartupMessages(library(org.Hs.eg.db))
                        genes_to_plot.first$Ensembl = mapIds(org.Hs.eg.db, 
                                                             keys = genes_to_plot.first$gene, 
                                                             keytype = "SYMBOL", 
                                                             column = "ENSEMBL")
                }
                
                # Parsing expression tabular data
                if(!identical(expression_tbl, F)){
                        exp_tbl = read.table(expression_tbl, header = T)
                        rna_seq_data = merge(x = genes_to_plot.first, 
                                             y = exp_tbl, 
                                             by.x = "Ensembl", 
                                             by.y = "fixed_ids", 
                                             all.x = T, 
                                             all.y = F)
                } else {
                        rna_seq_data = genes_to_plot
                }
                
                if(inv == TRUE){
                        # rotate/flip/invert genes within the inversion
                        
                        distance_from_start_st = (rna_seq_data[(rna_seq_data$start > inv_start 
                                                                & rna_seq_data$stop < inv_end),]$start
                                                  - inv_start)
                        
                        distance_from_start_en = (rna_seq_data[(rna_seq_data$start > inv_start 
                                                                & rna_seq_data$stop < inv_end),]$stop
                                                  - inv_start)
                        
                        
                        new_ends = (distance_from_start_en - inv_end )*-1
                        new_starts = (distance_from_start_st - inv_end )*-1
                        
                        rna_seq_data[(rna_seq_data$start > inv_start
                                      & rna_seq_data$stop < inv_end),"start"] = new_starts
                        
                        rna_seq_data[(rna_seq_data$start > inv_start
                                      & rna_seq_data$stop < inv_end),"stop"] = new_ends
                        
                }
                rna_seq_data$orientation = 1
                
                rna_seq_data$orientation[rna_seq_data$strand=="-"] = -1
                #### (b) plotting
                ########################################################################
                # First, if any expression data was provided
                if(!identical(expression_tbl, F)){
                        RNASeq_plot = 
                        ggplot(rna_seq_data, aes(xmin = start, xmax = stop, y = chrom, 
                                                 fill = log(TPM,2), color = log(TPM,2), forward = orientation)) +
                        geom_gene_arrow(arrowhead_height = unit(8, "mm"),
                                        arrowhead_width = unit(1,"mm"),
                                        arrow_body_height =  unit(8, "mm")) +
                        geom_vline(xintercept = del_start,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = del_end,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = inv_start,
                                   linetype = "dashed",
                                   color = "#b09c85") +
                        ggthemes::theme_few()+ scale_x_continuous(limits = x_range,
                                                                  expand = c(0, 0),
                                                                  labels=function(x)x/1e6) +
                        theme(axis.ticks.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.line.y.right = element_line(color = "black"),
                              axis.line.y.left = element_line(color = "black"),
                              axis.line.x.bottom = element_line(color = "black")) +
                        scale_y_discrete(expand = c(0,0.1)) +
                        scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA)+
                        scale_color_distiller(palette = "Blues", direction = 1, na.value = "gray80")+
                        guides(fill = "none", 
                               colour = guide_colourbar(direction = "vertical",
                                                        title = expression(log[2]~(TPM)))) +
                        xlab(paste0(chr_of_your_genes, " (Mb)")) +
                        ylab("Genes")
                } else {
                        # now, if no expression (expression_tbl FALSE)
                        RNASeq_plot = 
                        ggplot(rna_seq_data, aes(xmin = start, xmax = stop, y = chrom, forward = orientation)) +
                        geom_gene_arrow(arrowhead_height = unit(8, "mm"),
                                        arrowhead_width = unit(1,"mm"),
                                        arrow_body_height =  unit(8, "mm"),
                                        fill = NA, color = "gray80") +
                        geom_vline(xintercept = del_start,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = del_end,
                                   linetype = "dashed",
                                   color = "#00a087") +
                        geom_vline(xintercept = inv_start,
                                   linetype = "dashed",
                                   color = "#b09c85") +
                        ggthemes::theme_few()+ scale_x_continuous(limits = x_range,
                                                                  expand = c(0, 0),
                                                                  labels=function(x)x/1e6) +
                        theme(axis.ticks.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.line.y.right = element_line(color = "black"),
                              axis.line.y.left = element_line(color = "black"),
                              axis.line.x.bottom = element_line(color = "black")) +
                        scale_y_discrete(expand = c(0,0.1)) +
                        xlab(paste0(chr_of_your_genes, " (Mb)")) +
                        ylab("Genes")
                }
                
                ########################################################################
                ## Combined plot -- final\
                message("Here")
                final_plot = domainogram_plot / damid_plot / RNASeq_plot + 
                        plot_layout(heights = c(1,2.5,0.5),
                                    guides = "collect") &
                        theme(plot.margin = unit(c(0,0,0,0), "cm"))
                
                return(final_plot)

#                                                _
#                                             .-(_)
#                                            / _/
#                                         .-'   \
#                                        /       '.
#                                      ,-~--~-~-~-~-,
#                                     {__.._...__..._}             ,888,
#                     ,888,          /\##"  6  6  "##/\          ,88' `88,
#                   ,88' '88,__     |(\`    (__)    `/)|     __,88'     `88
#                  ,88'   .8(_ \_____\_    '----'    _/_____/ _)8.       8'
#                  88    (___)\ \      '-.__    __.-'      / /(___)
#                  88    (___)88 |          '--'          | 88(___)
#                  8'      (__)88,___/                \___,88(__)
#                            __`88,_/__________________\_,88`__
#                           /    `88,       |88|       ,88'    \
#                          /        `88,    |88|    ,88'        \
#                         /____________`88,_\88/_,88`____________\
#                        /88888888888888888;8888;88888888888888888\
#                       /^^^^^^^^^^^^^^^^^^`/88\\^^^^^^^^^^^^^^^^^^\
#                      /                    |88| \============,     \
#                     /_  __  __  __   _ __ |88|_|^  MERRY    | _ ___\
#                     |;:.                  |88| | CHRISTMAS, |      |
#                     |;;:.                 |88| |^  LISE     |      |
#                     |;;:.                 |88| '============'      |
#                     |::.                  |88|                     |
#                     |;;:'                 |88|                     |
#                     |:;,                  |88|             2022    |
#                     '---------------------""""---------------------'
        }
        
        
        
        
        #######################################################################
        ## Plot Domainograms starts here ######################################
        #######################################################################

        # changing the parameter names to be equal the previous script
        SMTH = smooth_window_size
        mygenes = strsplit(my_genes, split=",")[[1]]



        # changing the deletion start coorinate
        del_start = (inv_end - del_start) + inv_start
        Del<-data.frame(seqnames=recombination_chr,
                        start=c(del_start, del_end),
                        end=c(del_start+1, del_end+1))
        Inv<-data.frame(seqnames=recombination_chr,
                        start=c(inv_start, inv_end),
                        end=c(inv_start+1, inv_end+1))
        # png images resolution and size
        reso <- 300
        length <- 3.25*reso/72
        
        # Parse gene annotation
        logger("Parsing gene annotation")
        Refseq.mm10<-read.delim(annotation, head=FALSE, sep="\t", as.is=TRUE)
        Refseq.mm10<-Refseq.mm10[,c(3,5,6,13,12,4)]
        names(Refseq.mm10)<-c("chrom", "start", "stop", "gene", "score","strand")
        if(organism=="mouse"){ Refseq.mm10<-Refseq.mm10[-grep("Gm[0-9]*", Refseq.mm10$gene),] 
        } else { Refseq.mm10<-Refseq.mm10[-grep("LOC[0-9]*", Refseq.mm10$gene),] }
        logger(paste("Parsed",
                     nrow(Refseq.mm10), "genes; ",
                     length(unique(Refseq.mm10$gene)),
                     "unique genes"
        ))
        logger("-----------------------------------------------------------")
        logger("Finding your genes (--my_genes parameter)")
        ggene<-c(); gchrom<-c(); gstart<-c(); gstop<-c()
        for(g in mygenes){
                wh<-which(Refseq.mm10$gene==g)
                if(length(wh)>0){
                        ggene[g]<-as.character(g)
                        gchrom[g]<-as.character(unique(Refseq.mm10$chrom[wh]))
                        gstart[g]<-min(c(Refseq.mm10$start[wh],Refseq.mm10$stop[wh]))
                        gstop[g]<-max(c(Refseq.mm10$start[wh],Refseq.mm10$stop[wh]))
                }
                else(logger(paste("WARNING: cannot find gene: ", g)))
        }
        gchrom<-as.character(gchrom)
        logger("-----------------------------------------------------------")
        logger("Reading input data")
        ### Define input datasets
        Exp.Dfiles = strsplit(exp_dam_files, split=",")[[1]]
        Exp.DLfiles = strsplit(exp_antibody_files, split=",")[[1]]
        Control.Dfiles = strsplit(ctrl_dam_files, split=",")[[1]]
        Control.DLfiles = strsplit(ctrl_antibody_files, split=",")[[1]]
        
        data.summ<-paste0("Control (n = ", length(Control.DLfiles),"/", 
                          length(Control.Dfiles), ")")
        data.summ[2]<-paste0("Experimental (n = ", length(Exp.DLfiles),"/",
                             length(Exp.Dfiles), ")")
        logger(data.summ[1])
        logger(data.summ[2])
        
        #######################################################################
        # Performs the analysis only if it has not done before
        # (I wont indent this checkpoint for a better code readability)
        CHECKPOINT = paste0(output_dir, "/domainograms.RData")
        if(!file.exists(CHECKPOINT)){

        logger("-----------------------------------------------------------")
        logger("Combining data")
        #read and combine data:
        DamID.Control<-read.combine.DamID.replicates(DamLam=Control.DLfiles, 
                                                     Dam=Control.Dfiles)
        DamID.Exp<-read.combine.DamID.replicates(DamLam=Exp.DLfiles, 
                                                 Dam=Exp.Dfiles)
        
        #set extreme values of individual GATC fragments to the genome-wide mean:
        mcols(DamID.Control)$Dam[mcols(DamID.Control)$Dam > 100*mean(mcols(DamID.Control)$Dam)] <- mean(mcols(DamID.Control)$Dam)
        mcols(DamID.Control)$DamLam[mcols(DamID.Control)$DamLam > 100*mean(mcols(DamID.Control)$DamLam)]<-mean(mcols(DamID.Control)$DamLam)
        mcols(DamID.Exp)$Dam[mcols(DamID.Exp)$Dam > 100*mean(mcols(DamID.Exp)$Dam)] <- mean(mcols(DamID.Exp)$Dam)
        mcols(DamID.Exp)$DamLam[mcols(DamID.Exp)$DamLam > 100*mean(mcols(DamID.Exp)$DamLam)]<-mean(mcols(DamID.Exp)$DamLam)
        
        ##### INVERT FIRST
        if(inv == TRUE){
                logger("-----------------------------------------------------------")
                logger("Parsing the inversion (--recombination inversion)")
                # Invert sequences for Exp and for CTRL (since here, the control is also inverted)
                DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Inv$start[1] & end(DamID.Exp)< Inv$end[2]]$Dam <- rev(DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Inv$start[1] & end(DamID.Exp)< Inv$end[2]]$Dam)
                DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Inv$start[1] & end(DamID.Exp)< Inv$end[2]]$DamLam <- rev(DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Inv$start[1] & end(DamID.Exp)< Inv$end[2]]$DamLam)
                
                # VF230808: For the control now, do smoothing -> invertion 
		#           (as the control is now the wildtype) 
                # DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$Dam <- rev(DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$Dam)
                # DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$DamLam <- rev(DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$DamLam)
        }
        if(del == TRUE){
                logger("-----------------------------------------------------------")
                logger("Parsing the deletion (--recombination deletion)")
                # Delete the values within the deleted region, only for Exp
                Del.gr <- makeGRangesFromDataFrame(Del)
                ovl = findOverlaps(DamID.Exp, Del.gr)
                DamID.Exp <- GRangesList(DamID.Exp[1:(queryHits(ovl)[1]-1)], 
                                         DamID.Exp[queryHits(ovl)[1]:queryHits(ovl)[2]], 
                                         DamID.Exp[(queryHits(ovl)[2]+1):length(DamID.Exp)])
                DamID.Exp[[2]] <- NULL
                DamID.Exp <- unlist(DamID.Exp)
                
                
                # Create GRanges of the length of the deleted region (with NAs)
                deletion.gr <- DamID.Exp[queryHits(ovl)[1]:queryHits(ovl)[2]]
                deletion.gr$Dam <- NULL
                deletion.gr$DamLam <- NULL
                deletion.gr$logratio <- as.numeric(NA)
        }
        # if allele == "129" does nothing
        
        logger("-----------------------------------------------------------")
        logger("Smoothing")
        #smooth:
        sDamID.Control<-smooth.DamID(DamID.Control, ws=SMTH, pseud=30, title="control") #slightly higher pseudocount because of data sparsity
        sDamID.Exp<-smooth.DamID(DamID.Exp, ws=SMTH, pseud=30, title="experiment")
        # Combine the smoothed with the deleted GRanges
        if(del == TRUE) {
                logger("Combine the smoothed with the deleted data (--del TRUE)")
                sDamID.Exp <- c(sDamID.Exp[1:(queryHits(ovl)[1]-1)], deletion.gr,
                                sDamID.Exp[(queryHits(ovl)[1]):length(sDamID.Exp)])
        }
        
        logger("-----------------------------------------------------------")
        logger("Inspecting data")
        notna<-!(is.na(mcols(sDamID.Exp)$logratio) | is.na(mcols(sDamID.Control)$logratio))
        s<-sample(which(notna), 50000)
        par(mfrow=c(1,1))
        
        logger("-----------------------------------------------------------")
        logger("Creating scatter plots")
        logger("1) Raw (uncorrected) data")
        png(paste0(output_dir, "/scatter_ctrl_vs_exp_raw.png"), 
            units="in",res=reso,
            height=length,
            width=length)
        tmp = dev.cur()
        pdf(paste0(output_dir, "/scatter_ctrl_vs_exp_raw.pdf"))
        dev.control("enable")
        #### create pdf scatter_control_vs_exp
        plot(mcols(sDamID.Control)$logratio[s], mcols(sDamID.Exp)$logratio[s], 
             pch=19, cex=0.5, col="#00000022")
        #highlight Exp gene:
        COL<-rainbow(length(mygenes)) #plotting colors
        for(g in 1:length(mygenes)){
                w<-which(seqnames(sDamID.Exp)==gchrom[g] & start(sDamID.Exp)>=gstart[g] & end(sDamID.Exp)<=gstop[g])
                points(mcols(sDamID.Control)$logratio[w], mcols(sDamID.Exp)$logratio[w], pch=19, cex=0.5, col=COL[g])
                legend("topleft", pch=19, legend=mygenes, col=COL)
        }
        LM<-lm(mcols(sDamID.Exp)$logratio ~ mcols(sDamID.Control)$logratio)
        abline(LM, lwd=3, col="green")
        abline(a=0, b=1, col="black", lty="dashed", lwd=3)
        legend("bottomright", legend=c("diagonal", "linear fit"), 
               col=c("black", "green"), lwd=3, lty=c("solid", "dashed"))
        dev.copy(which=tmp)
        dev.off()
        dev.off()
        ### pdf ends here
        
        
        logger("Correlation coefficients for raw data:")
        logger(paste("Pearson:",
                     cor(mcols(sDamID.Control)$logratio, mcols(sDamID.Exp)$logratio, 
                         use="p",
                         method="p"),
                     "\nSpearman:",
                     cor(mcols(sDamID.Control)$logratio, mcols(sDamID.Exp)$logratio, 
                         use="p", 
                         method="s")
        ))
        #is slope close to 1:
        #logger(LM$coefficients)
        
        logger("Correcting with loess smoothing")
        LOESS50<-loess(mcols(sDamID.Exp)$logratio[s]-mcols(sDamID.Control)$logratio[s] ~ mcols(sDamID.Control)$logratio[s], 
                       span=0.5)
        #apply correction to experimental data:
        mcols(sDamID.Exp)$logratio <- mcols(sDamID.Exp)$logratio - predict(LOESS50, mcols(sDamID.Control)$logratio)
        #plot again: 
        logger("2) Corrected data")
        png(paste0(output_dir, "/scatter_ctrl_vs_exp_corrected.png"),
            units="in",res=reso,height=length,width=length)
        tmp1 = dev.cur()
        pdf(paste0(output_dir, "/scatter_ctrl_vs_exp_corrected.pdf"))
        dev.control("enable")
        ##### pdf 2 : after correction
        plot(mcols(sDamID.Control)$logratio[s], mcols(sDamID.Exp)$logratio[s], 
             pch=19, cex=0.5, col="#00000022")
        for(g in 1:length(mygenes))
        {w<-which(seqnames(sDamID.Exp)==gchrom[g] & start(sDamID.Exp)>=gstart[g] & end(sDamID.Exp)<=gstop[g])
        points(mcols(sDamID.Control)$logratio[w], mcols(sDamID.Exp)$logratio[w], 
               pch=19, cex=0.5, col=COL[g])
        legend("topleft", pch=19, legend=mygenes, col=COL)
        }
        LMcor<-lm(mcols(sDamID.Exp)$logratio[s] ~ mcols(sDamID.Control)$logratio[s])
        abline(LMcor, lwd=3, col="green")
        abline(a=0, b=1, col="black", lty="dashed", lwd=3)
        legend("bottomright", legend=c("diagonal", "linear fit"), 
               col=c("black", "green"), lwd=3, lty=c("solid", "dashed"))
        dev.copy(which=tmp1)
        dev.off()
        dev.off()
        #logger("after loess correction:")
        #logger(LMcor)
        
        logger("Correlation coefficients after correction:")
        logger(paste("Pearson:",
                     cor(mcols(sDamID.Control)$logratio, mcols(sDamID.Exp)$logratio,
                         use="p", method="p"),
                     "\nSpearman:",
                     cor(mcols(sDamID.Control)$logratio, mcols(sDamID.Exp)$logratio,
                         use="p", method="s")
        ))

        logger("-----------------------------------------------------------")
        logger("Calculating domainograms")

        #### INVERT FIRSTTTTT
        if(inv == TRUE){
                # Invert Control sequence
                # ^^^ is not necessary anymore, since the control here is already inverted
		# VF230808: now it is, because the control is the wildytpe

                DamID.Control[(seqnames(DamID.Control) == recombination_chr 
                               & start(DamID.Control) > Inv$start[1] 
                               & end(DamID.Control)< Inv$end[2])]$Dam <- 
                        rev(DamID.Control[(seqnames(DamID.Control) == recombination_chr 
                                           & start(DamID.Control) > Inv$start[1] 
                                           & end(DamID.Control)< Inv$end[2])]$Dam)
                DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$DamLam <- rev(DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Inv$start[1] & end(DamID.Control)< Inv$end[2]]$DamLam)
                
                
        }

        if(del == TRUE){
                DamID.Control<-read.combine.DamID.replicates(DamLam=Control.DLfiles, 
                                                             Dam=Control.Dfiles)
                DamID.Exp<-read.combine.DamID.replicates(DamLam=Exp.DLfiles, 
                                                         Dam=Exp.Dfiles)
                
                #set extreme values of individual GATC fragments to the genome-wide mean:
                mcols(DamID.Control)$Dam[mcols(DamID.Control)$Dam > 100*mean(mcols(DamID.Control)$Dam)] <- mean(mcols(DamID.Control)$Dam)
                mcols(DamID.Control)$DamLam[mcols(DamID.Control)$DamLam > 100*mean(mcols(DamID.Control)$DamLam)]<-mean(mcols(DamID.Control)$DamLam)
                mcols(DamID.Exp)$Dam[mcols(DamID.Exp)$Dam > 100*mean(mcols(DamID.Exp)$Dam)] <- mean(mcols(DamID.Exp)$Dam)
                mcols(DamID.Exp)$DamLam[mcols(DamID.Exp)$DamLam > 100*mean(mcols(DamID.Exp)$DamLam)]<-mean(mcols(DamID.Exp)$DamLam)
        }
        
        #calculate empirical domainograms:
        #w<-round(smooth_window_size*1.15^c(0:27)) %/%2 *2 +1 #odd windowsizes, exponentially increasing from ~40kb to ~2Mb
        # Changed here.
        # Now the first window of the domainogram corresponds to the DamID pair window (301 by default)
        # With this we would have less y-ticks in the domainograms
        # pos is the position in w of the first number of w bigger than your SMTH
        w = (round(67*1.15^c(0:27))%/%2 *2 +1)
        # minus 1 because we want to replace the previous position with your SMTH
        pos = Position(function(x) x > smooth_window_size, w) - 1 
        w = (round(67*1.15^c(0:27))%/%2 *2 +1)[pos:28]
        w[1] = smooth_window_size
        DG.Exp<-emp.domainogram(damid.exp=DamID.Exp, 
                                damid.ctl=DamID.Control, 
                                chrom=gchrom[1], 
                                window.sizes=w,
                                up=FALSE, 
                                logratios=TRUE, 
                                #loess.corr = NULL)
                                loess.corr=LOESS50)
        #DG.chr2<-emp.domainogram(damid.exp=DamID.Exp, damid.ctl=DamID.Control, chrom="chr2", window.sizes=w, up=FALSE, logratios=TRUE, loess.corr=LOESS50)
        
        if(TRUE) {
                # Clean the domainogram around the borders
                Del.gr <- makeGRangesFromDataFrame(Del)
                ovl.dom <- findOverlaps(DG.Exp, Del.gr)
                Inv.gr <- makeGRangesFromDataFrame(Inv)
                ovl.dom.inv <- findOverlaps(DG.Exp, Inv.gr)
        
                for(i in colnames(mcols(DG.Exp))){
                        g = as.integer(i)
                         
                        # INVERT FIRST
                        mcols(DG.Exp[(queryHits(ovl.dom.inv)[1]-(g+1)/2):(queryHits(ovl.dom.inv)[1]+(g+1)/2)])[i] <- 2
                        mcols(DG.Exp[(queryHits(ovl.dom.inv)[2]-(g+1)/2):(queryHits(ovl.dom.inv)[2]+(g+1)/2)])[i] <- 2

                        mcols(DG.Exp[(queryHits(ovl.dom)[1]-(g+1)/2):(queryHits(ovl.dom)[1]+(g+1)/2)])[i] <- 2
                        mcols(DG.Exp[(queryHits(ovl.dom)[2]-(g+1)/2):(queryHits(ovl.dom)[2]+(g+1)/2)])[i] <- 2
                }
        }
        

        ## INVERT FIRST
        if(inv == TRUE){
                # Invert smoothed sequences for Control only
                sDamID.Control[seqnames(sDamID.Control) == recombination_chr & start(sDamID.Control) > Inv$start[1] & end(sDamID.Control)< Inv$end[2]]$logratio <- rev(sDamID.Control[seqnames(sDamID.Control) == recombination_chr & start(sDamID.Control) > Inv$start[1] & end(sDamID.Control)< Inv$end[2]]$logratio)
                
                # NONE OF THIS IS PERFORMED HERE BECAUSE HERE THE CONTROL IS INVERTED
                # -------------------------------------------------------------------
                # VF230808: now it is, because the control is the wildytpe

                # remove values at the junction, only for Ctrl. To show that is 
                # the WT the lines are not connected.
                ovl = findOverlaps(sDamID.Control, Inv.gr)
                
                sDamID.Control[(queryHits(ovl)[1]-15):(queryHits(ovl)[1])]$logratio <- NA
                sDamID.Control[(queryHits(ovl)[1]):(queryHits(ovl)[1]+15)]$logratio <- NA
                
                sDamID.Control[(queryHits(ovl)[2]-15):(queryHits(ovl)[2])]$logratio <- NA
                sDamID.Control[(queryHits(ovl)[2]):(queryHits(ovl)[2]+15)]$logratio <- NA
        }

        if(del == TRUE){
                # Clean the domainogram within the deleted region
                mcols(DG.Exp[(queryHits(ovl.dom)[1]):(queryHits(ovl.dom)[2])]) <- 2
                
                logger("Number of read coming from deleted region:")
                logger(paste0("Experimental: ", 
                              sum(DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Del$start[1] & end(DamID.Exp)< Del$end[2]]$Dam, DamID.Exp[seqnames(DamID.Exp) == recombination_chr & start(DamID.Exp) > Del$start[1] & end(DamID.Exp)< Del$end[2]]$DamLam)
                ))
                # In Ctrl, non deleted
                logger(paste0("Control: ", 
                              sum(DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Del$start[1] & end(DamID.Control)< Del$end[2]]$Dam, DamID.Control[seqnames(DamID.Control) == recombination_chr & start(DamID.Control) > Del$start[1] & end(DamID.Control)< Del$end[2]]$DamLam)
                ))
        }
        
        

        # Save checkpoint
        logger(paste0("Saving a checkpoint for this analysis at: ", CHECKPOINT))
        save(sDamID.Exp, sDamID.Control, SMTH, data.summ,
             DG.Exp, mygenes, Refseq.mm10, Del, Inv, file = CHECKPOINT)
        
        } else {
                # If the analysis has already been done:
                logger("-----------------------------------------------------------", fg="green")
                logger(paste0("Checkpoint found: ", CHECKPOINT), fg="green")
                logger("Skipping all the calculations and jupping to plotting", fg="green")
                load(CHECKPOINT)}
        logger("-----------------------------------------------------------", fg="green")
        logger("Plotting domainograms")
        
        # Setting values for the parameters, based on the recombination
        h = 2100
        w = 3600
        reso = 300
        
        L_margin = 
                sprintf("%.1e", M_L) %>% 
                str_replace_all(c("\\+"="", "\\.0"="", "e0"="e"))
        
        R_margin = 
                sprintf("%.1e", M_R) %>% 
                str_replace_all(c("\\+"="", "\\.0"="", "e0"="e"))
        
        logger(paste0("Using left margin: ", L_margin, ", and right margin: ", R_margin))
        
         png(paste0(output_dir, 
                    "/domainogram_",
                    L_margin,
                    "L_", 
                    R_margin, 
                    "R.png"),
             res=reso,height=h,width=w)
         tmp2 = dev.cur()
         pdf(paste0(output_dir, 
                    "/domainogram_",
                    L_margin,
                    "L_", 
                    R_margin,
                    "R.pdf"), 
             width=12,height=7)
        dev.control("enable")
        p = plot_domainograms_with_ggplot(dg_exp = DG.Exp,
                                          Del = Del,
                                          Inv = Inv, 
                                          annot = Refseq.mm10, 
                                          my_genes = mygenes,
                                          recombination_centered = recombination_centered,
                                          mrg.left = M_L,
                                          mrg.right = M_R,
                                          smooth_window_size = smooth_window_size,
                                          sDamID.Exp = sDamID.Exp,
                                          sDamID.Control = sDamID.Control,
                                          logratio = T,
                                          y_limit = NA,
                                          expression_tbl = expression_tbl,
                                          organism = organism,
                                          plot_title = plot_title,
                                          del = del,
                                          inv = inv,
                                          del_start = del_start,
                                          del_end = del_end,
                                          inv_start = inv_start,
                                          inv_end = inv_end)
        #Sys.sleep(5)
        plot(p)
        #ggsave(paste0(output_dir,"/domainogram_", L_margin,"L_", R_margin,"R.pdf"),
        #       device = "pdf",
        #       width = 12,
        #       height = 7)
        dev.copy(which=tmp2)
        garbage = dev.off()
        garbage = dev.off()
        logger("-----------------------------------------------------------")
        logger("All done!", fg="green")
        logger(paste0("Domainogram can be found at: ", output_dir, 
                    "/domainogram_",
                    L_margin,
                    "L_", 
                    R_margin,
                    "R.pdf"), fg="green")
        logger("Bye", fg="green")
        
        
        
}



##############################################################################
## User input ################################################################
##############################################################################


# Get input
option_list <- list(
        make_option(c("--recombination"), type="character", default=NULL, 
                    help="Type of the recombination in this library (\'deletion\' or \'inversion\').\n\t\tIf it's a non-recombinant library, this should be \'FALSE\' (required)", 
                    metavar="inversion/deletion/FALSE"),
        make_option(c("--exp_dam_files"), type="character", default=NULL, 
                    help="Comma separated list of experimental Dam files in .txt.gz format (required)", 
                    metavar="list"),
        make_option(c("--exp_antibody_files"), type="character", default=NULL, 
                    help="Comma separated list of experimental antibody files in .txt.gz format (required)", 
                    metavar="list"),
        make_option(c("--ctrl_dam_files"), type="character", default=NULL, 
                    help="Comma separated list of control Dam files in .txt.gz format (required)",
                    metavar="list"),
        make_option(c("--ctrl_antibody_files"), type="character", default=NULL, 
                    help="Comma separated list of control antibody files in .txt.gz format (required)", 
                    metavar="list"),
        make_option(c("--annotation"), type="character", default=NULL, 
                    help="Annotation file (required)", metavar="path"),
        make_option(c("--output_dir"), type="character", default=NULL, 
                    help="Output directory -- will be created if not exists (required)", 
                    metavar="path"),
        make_option(c("--expression_tbl"), type="character", default=NULL, 
                    help="Path to expression data (check get_expression_tbl.R)", 
                    metavar="path"),
        
        # Optional arguments
        make_option(c("--organism"), type="character", default="mouse", 
                    help="Working organism (should be 'mouse' or 'human')", 
                    metavar="mouse/human"),
        make_option(c("--smooth_window_size"), type="integer", default=301, 
                    help="Smoothing window size (GATC fragments), must be an odd number (default: 301)", 
                    metavar="Int(odd)"),
        make_option(c("--my_genes"), type="character", default="Tdgf1", 
                    help="Comma separated list of genes to be the center of the plots (default: Tdgf1)", 
                    metavar="gene_name"),
        make_option(c("--inv_start"), type="integer", default=110927601, 
                    help="Start coordinate of the recombined region (default: 110927601)", 
                    metavar="coordenate"),
        make_option(c("--inv_end"), type="integer", default=112440815, 
                    help="End coordinate of the recombined region (default: 0112440815)", 
                    metavar="coordenate"),
                    make_option(c("--del_start"), type="integer", default=110927601, 
                    help="Start coordinate of the recombined region (default: 110927601)", 
                    metavar="coordenate"),
        make_option(c("--del_end"), type="integer", default=112440815, 
                    help="End coordinate of the recombined region (default: 0112440815)", 
                    metavar="coordenate"),
        make_option(c("--inv"), type="character", default="FALSE", 
                    help="Is there an inversion? (TRUE or FALSE)", 
                    metavar="true/false"),
        make_option(c("--del"), type="character", default="FALSE", 
                    help="Is there a deletion? (TRUE or FALSE)", 
                    metavar="true/false"),
        make_option(c("--recombination_chr"), type="character", default="chr9", 
                    help="Name of the recombined chromosome (default: chr9)", 
                    metavar="chromosome_name"),
        make_option(c("--recombination_centered"), type="character", default="FALSE", 
                    help="Should the plot be centered at the recombination? (TRUE or FALSE)\n\t\tIf FALSE, the center will be the genes in --my_genes", 
                    metavar="true/false"),
        make_option(c("--right_margin"), type="character", default=3e6, 
                    help="margin of the domainogram plot (default: 3e6)", 
                    metavar="int"),
        make_option(c("--left_margin"), type="character", default=3e6, 
                    help="Left margin of the domainogram plot (default: 3e6)", 
                    metavar="int"),
        make_option(c("--plot_title"), type="character", default="Domainogram", 
                    help="Title for the domainogram (default: Domainogram)", 
                    metavar="character")
)

opt_parser = OptionParser(usage = "\tRscript %prog [options]",
                          formatter = TitledHelpFormatter,
                          option_list=option_list,
                          description = paste0(
                                  "\n-----------------------------------------------------------",
                                  "\nPlot domainograms wrapper script",
                                  "\n-----------------------------------------------------------\n",
                                  "Version: ",
                                  script_version,
                                  "\nVinicius Franceschini-Santos, Lise Dauban 221207"),
                          epilogue = "Example of command line:\n------------------------\n\nRscript plot_domainograms.R \\\
       \t--exp_dam_files T48h20r16_1_Dam_CAS-gatc.counts.txt.gz \\\
       \t--exp_antibody_files T48h20r16_1_LB1_CAS-gatc.counts.txt.gz \\\
       \t--ctrl_dam_files T201_1_Dam_CAS-gatc.counts.txt.gz \\\
       \t--ctrl_antibody_files T201_1_LB1_CAS-gatc.counts.txt.gz \\\
       \t--annotation /DATA/usr/b.v.steensel/Projects/LADrewiring/data/ucsc/mm10/ncbiRefSeq.txt \\\
       \t--output_dir ./domainograms_output \\\
       \t--plot_title \"Clone T48h20r16\" \\\
       \t--recombination deletion \\\
       \t--organism mouse \\\
       \t--recombination_start 110927601 \\\
       \t--recombination_end 112877498 \\\
       \t--recombination_centered TRUE \\\
       \t--expression_tbl \"~/mydata/GitHub/plot_domainograms/4DNFIPFKK5LM_RNAseq_fixedIDs.tsv\"\n\n")
opt = parse_args(opt_parser)

# Check the allele parameter
# if(opt$allele != "129" & opt$allele != "CAST") {
#         print_help(opt_parser)
#         stop("The --allele parameter must be either \'129\' or \'CAST\'",call. = F)
# }

# Check the organism parameter
if(opt$organism != "mouse" & opt$organism != "human") {
        print_help(opt_parser)
        stop("The --organism parameter must be either \'mouse\' or \'human\'",call. = F)
}

# Check the recombination parameter for CAST
if((opt$recombination != "inversion" 
    & opt$recombination != "deletion"
    & opt$recombination != "FALSE")) {
        print_help(opt_parser)
        stop("The --recombination parameter must be either \'deletion\', \'inversion\', or \'FALSE\'",call. = F)
}

# Check required params
for(param in c(opt$exp_dam_files, opt$exp_antibody_files, opt$ctrl_dam_files, 
               opt$ctrl_antibody_files, opt$annotation, opt$output_dir,
               opt$recombination, opt$expression_tbl)){
        if(is.null(param)){
                stop("Not all required arguments present.\n", call.=FALSE)
        }
}

# Parse expression_tbl -- must be boolean if FALSE
if(opt$expression_tbl == "FALSE"){ 
        expression_tbl = as.logical(opt$expression_tbl) 
} else {
        expression_tbl = opt$expression_tbl
}


suppressMessages({
main(smooth_window_size = opt$smooth_window_size,
     my_genes = opt$my_genes,
     inv_start = opt$inv_start,
     inv_end = opt$inv_end,
     del_start = opt$del_start,
     del_end = opt$del_end,
     annotation = opt$annotation,
     exp_dam_files = opt$exp_dam_files,
     exp_antibody_files = opt$exp_antibody_files,
     ctrl_dam_files = opt$ctrl_dam_files, 
     ctrl_antibody_files = opt$ctrl_antibody_files,
     recombination_chr = opt$recombination_chr,
     output_dir = opt$output_dir,
     plot_title = opt$plot_title,
     script_version = script_version,
     inv = as.logical(opt$inv),
     del = as.logical(opt$del),
     recombination_centered = as.logical(opt$recombination_centered),
     organism = opt$organism,
     expression_tbl = expression_tbl,
     M_L = as.numeric(opt$left_margin),
     M_R = as.numeric(opt$right_margin))
})