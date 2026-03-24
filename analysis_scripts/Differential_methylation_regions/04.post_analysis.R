library(data.table)
library(dplyr)

mctypes <- c('BS','OXBS')


##################
for (mctype in mctypes){
    my_directory <- paste(mctype,"_DMR_results_ExSCZ_NoNorm/combined", sep = '')

    # Get a list of all .csv files in the directory
    csv_files <- list.files(path = my_directory, pattern = "\\.csv$")

    # Initialize an empty list to hold data frames
    list_of_dfs <- list()

    # Loop through each .csv file
    for (file_name in csv_files) {
    # Generate the full path to the .csv file
    full_path <- file.path(my_directory, file_name)
    
    # Read the .csv file into a data frame
    cat(full_path, '\n')

    df <- fread(full_path, sep = '\t')
    if(nrow(df) <= 1){

        next
    }
    df$chr_factor <- factor(sub("chr", "", df$chr), levels = as.character(c(1:22, "X")))

    sorted_df <- df %>% arrange(chr_factor, start)

    sorted_df$Effect <- ifelse(sorted_df$areaStat > 0,'-','+')
    sorted_bed <- sorted_df %>% select(chr, start,end,Effect)
    #sorted_bed <- sorted_bed[sorted_bed$chr != 'chrX']
    print(nrow(sorted_bed))
    output_name <- sub('csv','bed',full_path)
    fwrite(sorted_bed,output_name, row.names = F,col.names = F, sep = '\t')
    }
}