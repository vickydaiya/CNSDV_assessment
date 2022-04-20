#install.packages("openxlsx")

library(data.table)
library(openxlsx)

tsv_dir <- paste0(getwd(),"/Azimuth TSVs/")
tsv_files <- list.files(tsv_dir)

data_list <- list()
for (tsv_file in tsv_files)
{
  data = read.delim(paste0(tsv_dir,tsv_file), row.names = 1)
  data <- data.table(data)
  data_summary <- data[,list(count = .N, percentage = (.N/nrow(data))*100),by='predicted.cluster']
  data_list <- append(data_list,list(data_summary))
}
names(data_list) = tsv_files
write.xlsx(data_list, file = "Azimuth_analysis.xlsx",headerStyle = createStyle(textDecoration = "Bold"), colWidths = "auto")