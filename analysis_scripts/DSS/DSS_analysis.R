# Merge differentially methylated sites (DMS) to form regions (DMRs)

suppressMessages(library("DSS"))
suppressMessages(library(data.table))
suppressMessages(library(bsseq))
suppressMessages(library(dplyr))


args <- commandArgs(trailingOnly = TRUE)
group1 <- args[1]
group2 <- args[2]
celltype <- args[3]
Assay <- args[4]
qc = fread('/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/metadata_DNA_QC.tsv.gz')

DSS_analysis <- function(group1, group2, celltype, chr){
    cat('Process: ',group1, group2, celltype, chr,'\n')

    ### Load DML fit result
    file_name <- paste('/cndd/hex002/PsychEncode/methylation/Aged_DSS/',Assay,'_results_ExSCZ_NoNorm/',group1,'.',group2, '.',celltype, '_chr',chr,'_DSS.RData', sep = '')
    load(file_name)

    ### DML Test

    # dmlfit.age_group = DMLtest.multiFactor(dmlfit, coef=2)
    # dmr_result <- callDMR(dmlfit.age_group, p.threshold=0.05)

    ### Compute adjusted result
    # mc_group1 <-  mean(subset(qc, (DNA_passQC == 'True') & (period == group1) & (assay == Assay) & (celltype == celltype))$mCGlevel)
    # mc_group2 <-  mean(subset(qc, (DNA_passQC == 'True') & (period == group2) & (assay == Assay) & (celltype == celltype))$mCGlevel)

    # dmlTest <- callDML(dmlTest.sm, delta=0.1, p.threshold=0.05)
    # dmlTest <- subset(dmlTest,mu1 > 0.05 | mu2 > 0.05)
    # dmr <- callDMR(dmlTest, p.threshold=0.05)
    

    # output_name <- paste('/cndd/hex002/PsychEncode/methylation/Aged_DSS/',Assay,'_DMR_results_ExSCZ/', group1,'.',group2, '.',celltype, '_chr',chr,'_DMR.csv',sep = '')
    # if (!is.null(dmr)){
    # dmr$diff.Methy <- dmr$diff.Methy - (mc_group1 - mc_group2)
    # fwrite(dmr, output_name, row.names = F, sep = '\t')
    # }

    # ### Compute un-adjusted result
    dmlTest <- DMLtest.multiFactor(dmlTest.sm,coef = 2)
    dmr <- callDMR(dmlTest,p.threshold=0.05)

    if (group1 > group2){
         dmr$areaStat <- -(dmr$areaStat)

    }
    output_name <- paste('/cndd/hex002/PsychEncode/methylation/Aged_DSS/',Assay,'_DMR_results_ExSCZ_NoNorm/', group1,'.',group2, '.',celltype, '_chr',chr,'_DMR_unadjusted.csv',sep = '')
    if (!is.null(dmr)){
    fwrite(dmr, output_name, row.names = F, sep = '\t')
    }
    cat('Finish process',group1, group2, celltype, chr,'\n')
}


for (chr in c(1:22,'X')){
            
                DSS_analysis(group1,group2, celltype, chr)
            
}
