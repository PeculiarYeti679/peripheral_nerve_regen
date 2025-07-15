#preprocess data

#download libraries
library(readr) #read in tab delimited txt files
library(purrr) #functional programming; reduce()
library(dplyr) #data manipulation; summarise 
library(tools) #file name manipulation 

#get all files in the untarred CEL folder with .txt.gz file type
files <- list.files(".", pattern = "*.txt.gz", full.names = TRUE)

#check for that amount of files in the folder (30) in this case
length(files)

#read the files in and move to dataframes
data_list <- lapply(files, function(f) {
  #this is the sample name
  sample_name <- file_path_sans_ext(basename(f))
  df <- read_tsv(
    f,
    #skip the first 9 lines that do no contain actual data 
    skip = 9,
    show_col_types = FALSE,
    col_types = cols_only(
      ProbeName = col_character(),
      gProcessedSignal = col_double()
    )
  )
  
  # collapse duplicate ProbeIDs by averaging their expression
  df_clean <- df %>%
    group_by(ProbeName) %>%
    summarise(!!sample_name := mean(gProcessedSignal, na.rm = TRUE), .groups = "drop")
  
  #rename the "ProbeName" to the "ProbeID" so that all samples can be consistent
  colnames(df_clean)[1] <- "ProbeID"
  return(df_clean)
})
print(data_list)

# merge all samples into one matrix
expression_matrix <- reduce(data_list, ~ inner_join(.x, .y, by = "ProbeID"))



# move ProbeID to rownames
rownames(expression_matrix) <- expression_matrix$ProbeID
expression_matrix <- expression_matrix[, -1]

#save the files to reuse later
saveRDS(expression_matrix, file = "expression_matrix.rds")

