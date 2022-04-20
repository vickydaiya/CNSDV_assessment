#install.packages('Seurat')
library(Seurat)

current_dir <- getwd()

files <- list.files(path = current_dir)
for (file in files)
{
  if (endsWith(file,"mex.tar.gz"))
  {
    untar(tarfile = file, exdir = 'data/')
  }
}

data_dir = paste0(current_dir,"/data/")

file.rename(paste0(data_dir,"/filtered_feature_bc_matrix","/features.tsv"),paste0(data_dir,"/filtered_feature_bc_matrix","/genes.tsv"))
file.rename(paste0(data_dir,"/GW18_motor","/enes.tsv"),paste0(data_dir,"/GW18_motor","/genes.tsv"))
file.rename(paste0(data_dir,"/GW19_M1_all","/enes.tsv"),paste0(data_dir,"/GW19_M1_all","/genes.tsv"))
file.rename(paste0(data_dir,"/GW19_M1_CP","/enes.tsv"),paste0(data_dir,"/GW19_M1_CP","/genes.tsv"))


datafiles = list.files(path = data_dir)
for (datafile in datafiles)
{
  data <- Read10X(data.dir = paste0("data/",datafile,"/"))
  seurat_object <- CreateSeuratObject(counts = data)
  if (!file.exists("RDS_files"))
  {
    dir.create("RDS_files")
  }
  saveRDS(seurat_object, paste0("RDS_files/",datafile,".rds"))
}
rm(data)
rm(seurat_object)

