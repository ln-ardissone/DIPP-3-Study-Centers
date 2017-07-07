######################################
##DIPP 3 CITIES ANALYSIS SOURCE FILE##
##Lexi Ardissone
######################################


##############################################################################################
#NOTES ON CODE COMMENT FORMAT
##############################################################################################

#This source file has been organized into sections, where section titles begin with '#'.
#All function names are in all CAPS and indicated with '##'.
#Additional function comments describing the function are indicated with '###'.
#Description of arguments are denoted with ###ARGUMENTS### following their respective functions.


##############################################################################################
#REQUIRED PACKAGES
##############################################################################################

require(phyloseq)
require(ggplot2)
require(lubridate)
require(plyr)
require(reshape2)
require(gridExtra)
require(geepack)
require(vegan)
require(gdata)


##############################################################################################
#PREPROCESSING FUNCTIONS
##############################################################################################

###These functions are used to plot histograms of reads counts for taxa and samples;
###this should facilitate preprocessing decisions (e.g. removing samples/taxa with low read counts)

##TAXA_MEAN_VAR

###Takes a phyloseq object and calculates the mean and variance of taxa abundances at a specified taxonomic RANK.

taxa_mean_var = function(physeq, RANK) {
  #RANK indicates what taxonomic level to calculate mean & var
  dat = psmelt(physeq)
  form = as.formula(paste('~',RANK))
  mu = ddply(dat, form, function(x) mean(x$Abundance))
  sig2 = ddply(dat, form, function(x) var(x$Abundance))
  MV = merge(mu, sig2, by=RANK)
  names(MV)[2] = c('Mean')
  names(MV)[3] = c('Var')
  
  return(MV)  
}

###ARGUMENTS###
###physeq - a phyloseq object
###RANK - indicates what taxonomic rank you want to calculate mean and variance for

#.............................................................................................#

##MV_PLOT

###Plots mean against variance on a log10 scale; 
###plots the line mean=variance as a reference to assess overdispersion of data

mv_plot = function(MV, TITLE) {
  ggplot(MV, aes(Mean, Var)) + 
    geom_point() +
    scale_y_log10() + 
    scale_x_log10() + 
    geom_abline(intercept=0, slope=1, lty='dashed', size=3, color='#666666') + 
    ggtitle(TITLE) +
    ylab('Variance') + 
    annotate('text', x=1, y=1, label='mean=Var', angle=30, vjust=2, hjust=0)
}

###ARGUMENTS###
###MV - the output of 'taxa_mean_var'; a data.frame where rows are taxa of a specified taxonomic rank and columns are mean and variance


##############################################################################################
#METADATA DESCRIPTIVES
##############################################################################################

#STANDARD_ERROR

#Calculates standard error
se <- function(x) {
  if (length(x) == 1) {
    0
  } else {
    sd(x, na.rm=T)/sqrt(length(x))
  }
}

##############################################################################################
#DIVERSITY ANALYSES
##############################################################################################


#GET_RAREFACTION

##This function takes the arguments:
#df - which is a data frame where rows are samples and columns are species
#STEP - the increments to be subsampled...start with a number ~1/25 the sample size
##e.g. if most samples have ~100,000 reads, step=4000
##it uses vegan's rarecurve function
#Returns a dataframe with 3 colums --> use to merge with metadata & customize plots

get_rarefaction = function(df, STEP) {
    rar = rarecurve(df, step=STEP, label=FALSE)
    rar.comb = data.frame()
    for (i in seq_along(rar)) {
        df2 = data.frame(cbind(Species = rar[[i]],
        SampleSize=attr(rar[[i]], "Subsample")))
        df2$ID = unique(names(attr(rar[[i]], "Subsample")))[2]
        rar.comb = rbind(rar.comb, df2)
    }
    return(rar.comb)
}


##############################################################################################
#DIFFERENTIAL ABUNDANCES
##############################################################################################

##PREP_MOLTEN_TABLE

prep_molten_table = function(physeq, RANK, ABUND, ABUND_BY) {
  glom = tax_glom(physeq, RANK)
  dat = psmelt(glom)
  #dat = psmelt(physeq)
  
  #select abundant taxa within a particular variable (e.g. subject, age bin, etc.)
  #Abundance is based on median
  form = as.formula(paste('~',ABUND_BY, '+', RANK))
  keep.taxa = ddply(dat, form, function(x) median(x$Abundance) > ABUND)
  keep.taxa = keep.taxa[keep.taxa$V1 == TRUE,]
  num.kept = length(unique(keep.taxa[,2]))
  message('These parameters keep ', num.kept, ' taxa at the ', RANK, ' level.')
  
  
  dat.kept = subset(dat, dat[,RANK] %in% keep.taxa[,RANK])
  
  dat.kept
}

#ARGUMENTS:
#physeq is a phyloseq object
#RANK is the taxonomic level to evaluate (Phylum, Class, Order, etc.)
#ABUND is the abundance threshold; varies with type of data; uses MEDIAN 
#e.g. 0.001 on proportion data is equivalent to 0.1% of reads
#ABUND_BY gives the option to retain OTUs of a specified ABUND by a given variable; 
#e.g. Abundance varies by age, so I use age_bin...take OTUs of a minimum abundance in any given month.
###EXAMPLE USAGE: dipp.molten = prep_molten_table(dipp.filt.g.rar, 'Genus', 20, 'age_bin')

#.............................................................................................#

##GENERATE_SUBSETS

#this function is to be applied to molten dataframes
#The WINDOW_SIZE specifies the window range centered on the WINDOW_STEP
#WINDOW_STEP is easily speicifed with seq()
#RANK (arg in quotes) specifies the taxanomic level and is needed in the down-sampling chunk
#VALUE specifies the variable to take the median of in creating the subset/downsampling observations from the same subject that lie within a given time window. E.g. 'Abundance'

generate_subsets = function(molten, AGE, WINDOW_STEP, WINDOW_SIZE, RANK, VALUE) {
    molten$value = molten[,VALUE]
    molten$age = molten[,AGE]
    subsets = data.frame()
    for (i in WINDOW_STEP) {
        d = WINDOW_SIZE/2.0
        x1 = i - d
        x2 = i + d
        molten.sub = molten[molten[,'age'] >= x1,]
        molten.sub = molten.sub[molten.sub[, 'age'] <= x2,]
        molten.sub$window <- i
        if (nrow(molten.sub) < 10) {next; }
        subsets <- rbind(subsets, molten.sub)
    }
    
    #Within each window, there may be multiple samples from the same subject
    #Down sample by taking the median to account for repeated measures
    form = as.formula(paste('~window +', RANK, '+ mask_id'))
    ds = ddply(subsets, form, function(x) {
        #extract sample-specific variables
        c(value=median(x$value),
        age = median(x$age),
        age_at_sampling=median(x$age_at_sampling),
        current_status = median(x$current_status),
        total_reads = median(x$total_reads),
        total_reads_postproc = median(x$total_reads_postproc)
        )
    })
    
    #NOTE: sample_ids are lost in this process
    
    #extracts unique subject data
    ds.uniq = molten[!duplicated(molten$mask_id),]
    sub.names = c('mask_id', setdiff(names(ds.uniq), names(ds)))
    
    #matches unique subject data extracted in previous step with subject ID in extended, down-sampled, molten table
    #can take a few minutes to run
    #variables to retain can be changed in the following line (function will break if these variables are not present)
    dss = merge(ds,
    ds.uniq[,sub.names],
    #'Progressor_type'
    all.x=T,
    by = 'mask_id')
    
    return(dss)
}

###EXAMPLE USAGE: dipp.subsets = generate_subsets(dipp.molten, seq(91.5,732, by=30.5), 60, 'Genus', 'Abundance')

#.............................................................................................#

#TEST_MW

#returns a data.frame of Mann-Whitney results adjusted p.values by BH methods (default)
    ##within each data.frame and grouped by variable specified by FACET argument
#DS is a melted data frame
#VALUE is the name of the variable to tst for differences; e.g. corresponding to taxa value or alpha value
#AGE is the name of the variable that specifies the age bins within to perform each test (e.g. window)
#RANK is the taxonomic rank to analyze or any variable to iterate the test on
#GROUP is the variable identifying the groups to compare using Mann-Whitney (e.g. seroconverted); assumes logical values are provided
#FACET is the variable of subgroups to partition the data and perform MW test separately (e.g. site)

###EXAMPLE USAGE: test_mw(dipp.alpha.agebins, 'value', 'variable', 'seroconverted', 'site')


test_mw = function(DS, VALUE, AGE, RANK, GROUP, FACET) {
    #reassign variable names
    DS$rank = DS[,RANK]
    DS$group = DS[,GROUP]
    DS$value = DS[,VALUE]
    if (!is.null(FACET)) {
        DS$facet = DS[,FACET]
    }
    
    #DS$facet = DS[,FACET]

    #specify formula for ddply --> which variables to stratify on to perform MW test?
    if (is.null(FACET)) {
        form = as.formula(paste('~', RANK, '+', AGE))
    } else {form = as.formula(paste('~', FACET, '+', RANK, '+', AGE))}

    form2 = as.formula(paste('~', RANK))
 

    mw.p = ddply(DS, form, function(x)
    c(MW.p_value = wilcox.test(value~group, data=x, exact=F)$p.value,
    median.age_at_sampling = median(x$age_at_sampling),
    median.abund.TRUE = median(x[x$group == TRUE,]$value),
    se.abund.TRUE = se(x[x$group == TRUE,]$value),
    nTRUE = nrow(x[x$group == TRUE,]),
    median.abund.FALSE = median(x[x$group == FALSE,]$value),
    se.abund.FALSE = se(x[x$group == FALSE,]$value),
    nFALSE = nrow(x[x$group == FALSE,])
    ))

    #relabel column name indicating p.value
    if (!is.null(FACET)) {
        names(mw.p)[4] = paste('p_value')
    } else {names(mw.p)[3] = paste('p_value')}
    #calculate adjusted p-value; performs adjusment independent of facet (e.g.)
    mw.p$p.adjust = ddply(mw.p, form, function(x) p.adjust(x$p_value, method=c('fdr')))$V1
    
    mw.p2 = mw.p[complete.cases(mw.p$p.adjust),]
    
    #create 'significance' column to used for plotting functions
    mw.p2$significance <- 'insignificant'
    if (min(mw.p2$p.adjust, na.rm=TRUE) < 0.05) {
        mw.p2[mw.p2$p.adjust < 0.05,]$significance <- '0.05'
    }
    if (min(mw.p2$p.adjust, na.rm=TRUE) < 0.01) {
        mw.p2[mw.p2$p.adjust < 0.01,]$significance <- '0.01'
    }
    if (min(mw.p2$p.adjust, na.rm=TRUE) < 0.001) {
        mw.p2[mw.p2$p.adjust < 0.001,]$significance <- '0.001'
    }
    mw.p2$significance = as.factor(mw.p2$significance)
    
    return(mw.p2)
}


#.............................................................................................#


#TEST_GEE

#The `test_gee` function takes 
#a dataframe (dat) (melted, from `prep_molten_table`), 
#the variable (VAR) to be tested (e.g. `seroconverted`, `testx5`, `status`, etc.), and 
#the value (VALUE) of the response variable (e.g. `Abundance` for taxa analyses, `value` for diversity analyses); 
#only accommodates one variable at this time and uses the formula `Abundance~seroconverted*age_at_sampling`. 
#This function is most easily implemented using `ddply` on a molten data frame. 
#It returns a data frame of test results (`coef`) for each taxa. 
#P-values are adjusted using the `fdr` method as implemented by the `p.adjust` function.

#Example: model_status.div.6 = ddply(dipp.div, ~variable, function(x) test_gee6(x, 'status_rank', 'value'))

#The output of `test_gee` can be further used to identify significant effects. Essentially, it runs a bunch of GEEs and outputs the coefficients (and an adjusted p-value).

#SINGLE VARIABLE MODELS
#-------------------------------------------------------------------
test_age_gee = function(dat, VALUE) {
  form = as.formula(paste(VALUE, '~ I(mo_at_sampling^2)'))
  m <- geeglm(form, id=mask_id, data=dat)
  mm = coef(summary(m))
  names(mm)[4] = paste('p_value')
  mm$p.adjust = p.adjust(mm$p_value, method=c('fdr'))
  mm$effect = row.names(coef(summary(m)))
  return(mm)
}

test_sero_gee = function(dat, VALUE) {
  form = as.formula(paste(VALUE, '~ I(age_month^2)*seroconverted'))
  m <- geeglm(form, id=mask_id, data=dat)
  mm = coef(summary(m))
  names(mm)[4] = paste('p_value')
  mm$p.adjust = p.adjust(mm$p_value, method=c('fdr'))
  mm$effect = row.names(coef(summary(m)))
  return(mm)
}

test_adj_gee = function(dat, VALUE) {
    form = as.formula(paste(VALUE, '~ I(age_month^2)*seroconverted*site + I(age_month^2)*Gender + receiving_breast_milk + antibiotics'))
    m <- geeglm(form, id=mask_id, data=dat)
    mm = coef(summary(m))
    names(mm)[4] = paste('p_value')
    mm$p.adjust = p.adjust(mm$p_value, method=c('fdr'))
    mm$effect = row.names(coef(summary(m)))
    return(mm)
}
#___________________


#.............................................................................................#
#GET_GEE_SIG_TAXA

#Extraacts significant results from test_XX_gee into a data frame

get_gee_sig_taxa = function(TEST_GEE, CUTOFF, EFFECT) {
    x = TEST_GEE
    y = as.numeric(CUTOFF)
    x.effect = subset(x, effect %in% EFFECT)
    #form2 = as.formula(paste(molten, '$', RANK))
    keep.sig = x.effect[x.effect$p.adjust <= y,]
    #sig.taxa = unique(keep.sig[,1])
    #taxa.prop.sig = subset(molten, form2 %in% sig.taxa)
    return(keep.sig)
}

#Example: modelB.div.6.sig = get_gee_sig_taxa(modelB.div.6, 0.05)

#.............................................................................................#
#PLOT_GEE_SIG_TAXA
plot_gee_sig_taxa = function(SIG, DS, VALUE, RANK, LTY, COLOR_VAR, FACET, TIMEPOINT) {
    #SIG is a vector of taxa names of interest --> allows you to plot 1 taxa at a time
    #DS is the melted data.frame including taxa abundances and matadata
    #VALUE is the name of the variable indicating taxa values (e.g. Abundance, value, etc.)
    #RANK is the taxa level to analyze
    #LTY is the the variable name to correspond to plot line type
    #COLOR_VAR is the the variable name to correspond to plot line color
    #FACET is the variable name to split the data by and make a separate plot for each (e.g. site)
    #TIMEPOINT is the age variable to use; name may vary depending on how samples were binned (e.g. age_month)
        #ASSUMES timepoint variable is specified
    
    
    #SETUP:
    #extract taxa of interest from melted data.frame
    
    taxa.df1 = DS
    taxa.df1$rank = as.factor(taxa.df1[,RANK])
    taxa.df = taxa.df1[taxa.df1$rank == SIG,]
    
    #relabel variables
    taxa.df$facet = as.factor(taxa.df[,FACET])
    taxa.df$lty = as.factor(taxa.df[,LTY])
    taxa.df$color_var = as.factor(taxa.df[,COLOR_VAR])
    taxa.df$value = taxa.df[,VALUE]
    taxa.df$timepoint = taxa.df[,TIMEPOINT]
    #assign discrete timepoints (monthly)
    #define formulas
    form = as.formula(paste('~', RANK))

    #Set up rug colors
    sig.colors = c(insignificant='gray')
    sig.colors <- c(sig.colors, '0.05'="blue")
    sig.colors <- c(sig.colors, '0.01'="green")
    sig.colors <- c(sig.colors, '0.001'="brown")
    sig.colors = c(sig.colors, 'TRUE'='#E84646')
    sig.colors = c(sig.colors, 'FALSE'='#65ADC2')
    
    sig.scale <- scale_colour_manual(values=sig.colors)
    
    #get significance values
    ds.MW = test_mw(taxa.df, VALUE, TIMEPOINT, RANK, COLOR_VAR, FACET)
    ds.MW$rank = ds.MW[,RANK]
    ds.MW$facet = ds.MW[,FACET]
    
    #Merge melted abundance table AND MW results
    diff_sig2 = merge(taxa.df, ds.MW, by=c('rank', TIMEPOINT, 'facet'), all.x=TRUE)

    #add age of sc
    meta = taxa.df[!duplicated(taxa.df$mask_id),]
    sc.ages = ddply(meta, ~facet, function(x) c(summary(x$age_first_sc)/30.5))
    names(sc.ages) = c('facet', 'Min', 'IQR1', 'Median', 'Mean', 'IQR3', 'Max', 'NA')
    diff_sig = merge(diff_sig2, sc.ages, by='facet', all.x=TRUE)

    #PLOT
    p = ggplot(diff_sig, aes(timepoint, value, lty=lty, color=as.logical(as.numeric(color_var)-1))) +
    geom_smooth(method='loess', alpha=0.2, size=2) +
    geom_vline(aes(xintercept=diff_sig$IQR1), lty='dashed') +
    geom_vline(aes(xintercept=diff_sig$IQR3), lty='dashed') +
    geom_vline(aes(xintercept=diff_sig$Median)) +
    facet_wrap(~facet) +
    scale_x_continuous(breaks=seq(3,24, by=3), labels=seq(3,24, by=3)) +
    geom_rug(aes(color=significance, size=1), sides='b') +
    theme_bw() +
    sig.scale +
    ylab('Relative abundance (proportion)') +
    xlab('Age at sampling (months)') +
    theme(legend.position='none')
    
    return(p)
}


plot_gee_sig_taxa_var = function(SIG, DS, VALUE, RANK, COLOR_VAR, TIMEPOINT, LINE_COLOR) {
    #SIG is a vector of taxa names of interest --> allows you to plot 1 taxa at a time
    #DS is the melted data.frame including taxa abundances and matadata
    #VALUE is the name of the variable indicating taxa values (e.g. Abundance, value, etc.)
    #RANK is the taxa level to analyze
    #COLOR_VAR is the the variable name to correspond to plot line color
    #TIMEPOINT is the age variable to use; name may vary depending on how samples were binned (e.g. age_month)
        #ASSUMES timepoint variable is specified
    
    
    #SETUP:
    #extract taxa of interest from melted data.frame
    
    taxa.df1 = DS
    taxa.df1$rank = as.factor(taxa.df1[,RANK])
    taxa.df = taxa.df1[taxa.df1$rank == SIG,]
    
    #relabel variables
    #taxa.df$facet = as.factor(taxa.df[,FACET])
    #taxa.df$lty = as.factor(taxa.df[,LTY])
    taxa.df$color_var = as.factor(taxa.df[,COLOR_VAR])
    taxa.df$value = taxa.df[,VALUE]
    taxa.df$timepoint = taxa.df[,TIMEPOINT]
    #assign discrete timepoints (monthly)
    #define formulas
    form = as.formula(paste('~', RANK))

    #Set up rug colors
    sig.colors = c(insignificant='gray')
    sig.colors <- c(sig.colors, '0.05'="blue")
    sig.colors <- c(sig.colors, '0.01'="green")
    sig.colors <- c(sig.colors, '0.001'="brown")
    sig.colors = c(sig.colors, 'TRUE'=LINE_COLOR[1])
    sig.colors = c(sig.colors, 'FALSE'=LINE_COLOR[2])
    
    sig.scale <- scale_colour_manual(values=sig.colors)

    #get significance values
    ds.MW = test_mw(DS, VALUE, TIMEPOINT, RANK, COLOR_VAR, NULL)
    ds.MW$rank = ds.MW[,RANK]
    #ds.MW$facet = ds.MW[,FACET]

    #Merge melted abundnace table AND MW results
    diff_sig2 = merge(taxa.df, ds.MW, by=c('rank', TIMEPOINT), all.x=TRUE)

    #add age of sc
    meta = taxa.df[!duplicated(taxa.df$mask_id),]
    sc.ages = ddply(meta, ~rank, function(x) c(summary(x$age_first_sc)/30.5))
    names(sc.ages) = c('rank', 'Min', 'IQR1', 'Median', 'Mean', 'IQR3', 'Max', 'NA')
    diff_sig = merge(diff_sig2, sc.ages, by='rank', all.x=TRUE)

    #PLOT
    p = ggplot(diff_sig, aes(timepoint, value, color=as.logical(as.numeric(color_var)-1))) +
    geom_smooth(method='loess', alpha=0.2, size=2) +
    geom_vline(aes(xintercept=diff_sig$IQR1), lty='dashed') +
    geom_vline(aes(xintercept=diff_sig$IQR3), lty='dashed') +
    geom_vline(aes(xintercept=diff_sig$Median)) +
    #facet_wrap(~facet) +
    scale_x_continuous(breaks=seq(3,21, by=3), labels=seq(3,21, by=3)) +
    geom_rug(aes(color=significance, size=1), sides='b') +
    theme_bw() +
    sig.scale +
    ylab('Relative Abundance (proportion)') +
    xlab('Age at sampling (months)')
    #theme(legend.position='none')
    
    return(p)
}

#.............................................................................................#

