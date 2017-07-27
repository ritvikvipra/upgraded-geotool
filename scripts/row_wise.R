
#==========Libraries required============
library(reshape)
library(reshape2)
library(tidyr)
library(tibble)
library(dplyr)

#==Define input file for data extraction==
wd <- "C:/Users/ritvik.vipra/Documents/misc"
setwd(wd)
input_filename <- "Tagged_GSM_Details_GSEfromMRV_dk02062017"

#==Load intermediate data into R==
data <- read.csv(paste0(input_filename, ".csv"), header = T, stringsAsFactors = F)
working_data <- data


#==Preprocess the data if needed (will change for each dataset)==
working_data$characteristics_ch1[22087:26408] <- sapply(working_data$characteristics_ch1[22087:26408],
                                                        FUN = function(x) paste0("strain: ", x))

#==Identify unique GSE's==
#unique_gse <- unique(working_data$GSE)

#==Store exact names in order of desired columns to be extracted from "characteristics_ch1" column==
final_des_col <- c("disease", "cell/cell line","cell type", "genotype/variationnn", "strainnn/background",
                   "tissue")
#==Create search libraries(update them as needed)==
desired_col_lib <- c("disease", "disease state", "disease status", "cell", "cell line", 
                     "cell type", "cell lineages", "genotype/variation","genotype", "strain/background", 
                     "strain", "background strain", "tissue", "Tissue", "tissue type", "tissue origin")
des_col_lib_disease <- c("disease", "disease state", "disease status","tumor type", "diagnosis", "disease.phenotype")
des_col_lib_cell <- c("cell", "cell line")
des_col_lib_celltype <- c("cell type")
des_col_lib_geno <- c("genotype/variation","genotype")
des_col_lib_strain <- c("strain/background", "strain", "background strain")
des_col_lib_tissue <- c("tissue", "Tissue", "tissue type", "tissue origin", "tissue_type", "organ/tissue")

#==List of above search libraries==
check_list <- list(des_col_lib_disease, des_col_lib_cell,des_col_lib_celltype, des_col_lib_geno,
                   des_col_lib_strain, des_col_lib_tissue) 

#==Create a dataframe to store appended results for all GSE id's==
#==Rows in this dataframe will be indexed(primary key) by (GSE id,GSM id) tuples==
m <- matrix(nrow = 0, ncol = 9)
final_df <- as.data.frame(m)
colnames(final_df) <- c("GSE","gsm", final_des_col, "Others")
rm(m)


#==Helper function definitions==
# function 2 - custom dataframe appending function
append_df <- function(main_df, small_df){
  main_df <- rbind(main_df, small_df)
  return(main_df)
}

# function - Logic for splliting string
create_aux_df <- function(charVec){
  charVecOdd <- as.vector(charVec[c(TRUE,FALSE)])
  charVecEven <- as.vector(charVec[c(FALSE,TRUE)])
  if(length(charVecOdd) > length(charVecEven)){
    charVecEven[length(charVecOdd)] <- "NA"
  }
  m <- matrix(nrow = 0, ncol = length(charVecOdd))
  tempDF <- as.data.frame(m)
  rm(m)
  tempDF <- rbind(tempDF, charVecEven)
  colnames(tempDF) <- charVecOdd
  return(tempDF)
}

# function 3 - For generating "Others" string from rejected columns
concatenate <- function(r){
  others_string <- "" 
  p <- paste(cols_to_others, as.character(r[cols_to_others]), sep = ":")
  for(k in p){
    if(nchar(others_string) == 0){others_string = paste(others_string, k, sep="")}
    else{others_string = paste(others_string, k, sep=";") }
  }
  if(others_string == ""){return(NA)}
  else{return(others_string)}
}

#===Columns for which to extract data initially=== 
col_names <- c("GSE", "gsm", "characteristics_ch1")

#===Actual data extraction module starts here===
start_time <- Sys.time()
for(row in seq(nrow(working_data)))
{
  # Getting the 'col_names' subsetted data for current row 
  curr_df <- subset(working_data[row, ], select = col_names)
  
  # Processing rows for which "characteristics_ch1" = NA
  if(is.na(as.character(curr_df$characteristics_ch1)))
  {
    na_row <- rep(NA, 9)
    final_df <- rbind(final_df, na_row)
    colnames(final_df) <- c("GSE","gsm", final_des_col, "Others")
    final_df[nrow(final_df),"GSE"] <- curr_df$GSE
    final_df[nrow(final_df),"gsm"] <- curr_df$gsm
    paste0("Completed for ", curr_df$GSE,"-",curr_df$gsm,". (all are NA)")
  }
  
  # Processing non-NA "characteristics_ch1" rows 
  if(!(is.na(as.character(curr_df$characteristics_ch1))))
  {
    #Preparing auxilliary df by splitting string
    tempList <- strsplit(curr_df$characteristics_ch1, split = "[;:]")
    DF <- as.data.frame(lapply(tempList, FUN = create_aux_df)[1])
    DF[] <- lapply(DF, as.character)
    aux_df <- DF
    aux_df <- cbind(curr_df$gsm, aux_df)
    aux_df <- cbind(curr_df$GSE, aux_df)
    colnames(aux_df)[1] <- "GSE"
    colnames(aux_df)[2] <- "gsm"
    aux_df[] <- lapply(aux_df, as.character)
    all_col <- colnames(aux_df)[-1][-1]
    Others <- rep("NA", times = nrow(aux_df))
    aux_df <- cbind(aux_df, Others)
    rm(Others)
    
    is_col_desired <- vector(mode = "numeric", length = length(all_col))
    
    # Main search and pattern matching algorithm
    iter <- 1
    for(k in check_list)
    {
      v = as.vector(sapply(all_col, FUN = function(i)
      {return(tryCatch((ifelse(any(grepl(i, k, ignore.case = T)), TRUE, FALSE)),
                       error = function(e) NA))}))
      if(any(v)) 
      {
        is_col_desired[which(v)] = 1
        all_col[which(v)] = final_des_col[[iter]]
      }
      iter <- iter + 1
    }
    
    # Determing cols to keep and cols to put into "Others"
    colnames(aux_df)[seq(3, ncol(aux_df)-1)] <- all_col 
    cols_to_take <- all_col[is_col_desired == 1]
    cols_to_others <- all_col[is_col_desired == 0]
    
    # Generating the "Others" string in original "characteristics_ch1" format
    complete_Others <- as.vector(apply(aux_df, 1, FUN = concatenate))
    aux_df$Others <- complete_Others
    aux_df <- aux_df[, !(names(aux_df) %in% cols_to_others)]
    n <- colnames(aux_df)[seq(3, length(colnames(aux_df))-1)]
    
    # Determining which cols in final_des_col are not present and adding them with NA value
    # in appropriate position and standardizing 'aux_df' to 'aux_df_mod' to add to 'final_df'
    f <- final_des_col %in% n
    final_des_col[f]
    append_list <- lapply(final_des_col[!(f)], FUN = function(x) x = rep(NA, nrow(aux_df)))
    names(append_list) <- final_des_col[!(f)]
    aux_df_mod = aux_df[,c("GSE","gsm", final_des_col[f], "Others")]
    iter <- 1
    for(l in as.tibble(append_list))
    {
      aux_df_mod <- add_column(aux_df_mod,l, .after = 2)
      colnames(aux_df_mod)[3] <- names(as.tibble(append_list))[iter]
      iter = iter + 1
    }
    aux_df_mod <- aux_df_mod[,c("GSE","gsm", final_des_col, "Others")]
    
    # Appending aux_df_mod to final_df
    final_df <- append_df(final_df, aux_df_mod) 
    print(paste0("Completed for ", curr_df$GSE,"-", curr_df$gsm , "."))
    
  }
}
end_time <- Sys.time()
total_runtime <- end_time - start_time
print(total_runtime)

#===Combine final_df with data===
output_data <- working_data
# selecting cols from final_df to add to output_data
final_cols_to_output <- final_df[,3:ncol(final_df)]
# Adding desired cols in proper positions with correct names
for(j in ncol(final_cols_to_output):1)
{
  output_data <- add_column(output_data, as.vector(final_cols_to_output[,j]), .after = 11)
  colnames(output_data)[12] <- colnames(final_cols_to_output)[j]
}

#===Set filepath of location where to store output===
output_path <- "C:/Users/ritvik.vipra/Documents"
write.csv(output_data, paste0(output_path, "/finaloutput.csv"))














