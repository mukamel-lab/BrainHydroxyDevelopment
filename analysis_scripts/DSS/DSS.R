suppressMessages(library(DSS))
source('modDSS.R')
suppressMessages(library(data.table))
suppressMessages(library(bsseq))
suppressMessages(library(dplyr))
####

args <- commandArgs(trailingOnly = TRUE)
Group1 <- args[1]
Group2 <- args[2]
Celltype <- args[3]
Assay <- args[4]


##
Sample2include = c('True') 
### Read Info

qc <- fread('/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/metadata_DNA_QC.tsv.gz')
data_dir <- '/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/allcg_rmSNPs/'

read_custom <- function(file) {
    df <- fread(file, select = c(1, 2, 5, 6))
    df <- df[, .(V1, V2, V6, V5)]
    
    # Add 'chr' prefix to the first column
    df$V1 <- ifelse(grepl("^[0-9XY]+$", df$V1), paste0("chr", df$V1), df$V1)
    
    setnames(df, c("chr", "pos", "N", "X"))
    return(df)
}

# Read files

###
group1_ids <- subset(qc, (DNA_passQC %in% Sample2include) & (period == Group1) & (assay == Assay) & (celltype == Celltype))$donor
group1_sex <- subset(qc, (DNA_passQC %in% Sample2include) & (period == Group1) & (assay == Assay) & (celltype == Celltype))$sex

# Create filenames using vectorized paste
filenames_1 <- paste0(data_dir, group1_ids, '_',Celltype,'_',Assay,'_hisat3n_rmSNPs_rmblacklist.CGN-Merge.allc.tsv.gz')

# Pre-allocate memory for data_list_1
num_files <- length(filenames_1)
data_list_1 <- vector("list", length = num_files)

# Read files with error handling
for (i in seq_along(filenames_1)) {
    data_list_1[[i]] <- tryCatch({
        read_custom(filenames_1[i])
    }, error = function(e) {
        cat("Error with Sample ID:", sample_ids[i], "\n")
        return(NULL)
    })
}

# Remove NULL entries
data_list_1 <- data_list_1[sapply(data_list_1, function(x) !is.null(x))]

# Efficiently name list elements
names(data_list_1) <- paste0("dat_1.", seq_along(data_list_1))


# Filtering data for Group2
group2_ids <- subset(qc, (DNA_passQC %in% Sample2include) & (period == Group2) & (assay == Assay) & (celltype == Celltype))$donor
group2_sex <- subset(qc, (DNA_passQC %in% Sample2include) & (period == Group2) & (assay == Assay) & (celltype == Celltype))$sex

# Create filenames using vectorized paste
filenames_2 <- paste0(data_dir, group2_ids, '_',Celltype,'_',Assay,'_hisat3n_rmSNPs_rmblacklist.CGN-Merge.allc.tsv.gz')

num_files <- length(filenames_2)
data_list_2 <- vector("list", length = num_files)

# Read files with error handling
for (i in seq_along(data_list_2)) {
    data_list_2[[i]] <- tryCatch({
        read_custom(filenames_2[i])
    }, error = function(e) {
        cat("Error with Sample ID:", sample_ids[i], "\n")
        return(NULL)
    })
}

# Remove NULL entries
data_list_2 <- data_list_2[sapply(data_list_2, function(x) !is.null(x))]

# Efficiently name list elements
names(data_list_2) <- paste0("dat_2.", seq_along(data_list_2))







############################################# Merge the data and metadata ############
Groups <- c(rep(Group1,length(data_list_1)),rep(Group2,length(data_list_2)))
Sex <- c(group1_sex,group2_sex)  
design = data.frame(Groups,Sex)
merged_list <- c(data_list_1, data_list_2)

######## RUN DSS per Chromsome
chromosomes <- c(paste0("chr", 1:22), "chrX")
for (chrom in chromosomes){

                # Split by chromosome
                chrom_sub_list <- lapply(merged_list, function(df) {
                return(df[df$chr == chrom])
                })

                BSobj <- makeBSseqData(chrom_sub_list,c(names(data_list_1),names(data_list_2)))  # Adjust the vector as needed
                print('all files loaded. Perform DSS:')

                ####  general experimental design
                dmlTest.sm =  DMLfit.multiFactor(BSobj, design=design, formula=~Groups+Sex, smoothing=TRUE)

                #### No covariates
                #dmlTest.sm = DMLtest(BSobj, group1=names(data_list_1), group2=names(data_list_2),smoothing=TRUE)

                # Junhao's Modfied Version
                # dmlTest.sm <- DMLtestMod(BSobj, smoothing = T,
                                #   group1 = names(data_list_1),
                                #   group2 = names(data_list_2))

                outputname <- paste('/cndd/hex002/PsychEncode/methylation/Aged_DSS/',Assay,'_results_ExSCZ_NoNorm/', Group1, '.', Group2,'.' ,Celltype,'_',chrom,'_DSS.RData', sep = '')
                save(dmlTest.sm, file = outputname)

                print(paste(chrom,'Done.'))
}

print('ALL Chromosome Done')


