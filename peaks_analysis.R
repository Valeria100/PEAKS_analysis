# 1 - Read samples file ----

#     Read the file containing all the info related to the samples: cell lines or tissue
samples_file <- read.csv(paste(getwd(),"samples_file.csv",sep="/"))
levels(samples_file$tissue)
#     Group the samples per tissue/cell type
samples_groups <- lapply(levels(samples_file$tissue), function(i) samples_file[which(samples_file$tissue==i),])
#     Read all the file names
file_names <- list.files(paste(getwd(),"/peaks",sep=""),full.names = TRUE)
#     Keep only the .narrowPeaks (already bedfile)
file_names_bed <- file_names[grep("_peaks.xls",file_names)]
#     Get the corresponding file path to the file name
samples_corr_filepath <- lapply(samples_groups, function(i) 
                                 sapply(i$ID, function(j) file_names_bed[grep(j,file_names_bed)]))
#     Substitute an empty cell with NA
samples_corr_filepath <- lapply(samples_corr_filepath, function(i)
                                sapply(i, function(j) ifelse(length(j)==0,NA,j)))
#     Merge the samples_group with the corresponding file path
samples_info <- lapply(1:length(samples_groups), function(i) cbind(samples_groups[[i]],samples_corr_filepath[[i]]))
#     Rename the columns for later
for(i in 1:length(samples_info)){
  colnames(samples_info[[i]]) <- c(colnames(samples_info[[1]])[-length(colnames(samples_info[[1]]))],"filePath")
}
#     Read files
samples_group_data <- vector("list", length(samples_info))
for(i in 1:length(samples_info)){
  samples_group_data[[i]] <- vector("list", nrow(samples_info[[i]]))
  for(j in 1:nrow(samples_info[[i]])){
    if(is.na(samples_info[[i]]$filePath[j])){
      samples_group_data[[i]][[j]] <- NA
    }else{
      samples_group_data[[i]][[j]] <- read.table(as.character(samples_info[[i]]$filePath[j]),skip=26,header=TRUE)
    }
  }
}

# 2 - Select the reliable peaks in each group ----

#     Select the peaks based on pvalue threshold (the value is -log10(pvalue)) 
# NOT FOR NOW

#     Intersect and union the peaks in each group using BEDTOOLS

write.table(samples_group_data[[1]][[1]][,c(1,2,3)],file="temp1.bed", quote=FALSE, col.names=FALSE, row.names=FALSE,sep="\t")
write.table(samples_group_data[[1]][[2]][,c(1,2,3)],file="temp2.bed", quote=FALSE, col.names=FALSE, row.names=FALSE,sep="\t")

system("bedtools intersect -a temp1.bed -b temp2.bed -wa -wb > intersected_temp12.bed")




intersect_different_windows <- function(MP, win){
  MP$start <- MP$start-win
  # MP$start[which(MP$start<1)] <- 0 
  MP$end <- MP$start+win
  
  MP <- cbind(MP,c(1:nrow(MP)))
  
  write.table(MP, file=paste("merged_peaks_window_",win,".bed",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
  
  system(paste("bedtools intersect -a merged_peaks_window_",win,
               ".bed -b tfs_location.bed -wa -wb > intersected_peaks_tfs_",win,".bed",
               sep=""))
  
  intersected_peaks <- read.table(paste("intersected_peaks_tfs_",win,".bed",sep=""), header=FALSE)
  intersected_peaks_noind <- get_minmax_range(intersected_peaks[,c(1,2,3,6,7)], "intersected_genes_noind.bed")
  intersected_peaks <- cbind(intersected_peaks_noind, intersected_peaks[,8])
  write.table(intersected_peaks, file=paste("intersected_tfs_",win,".bed",sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
  
  length(unique(intersected_peaks[,4]))
  
  ind_ar_tfs <- sort(unique(intersected_peaks[,4]))
  
  #TFs in around the AR peaks - directly regulated
  tfs_ar <- ordered_tfs[ind_ar_tfs,]
  return(tfs_ar)
}


tfs_ar_1500 <- intersect_different_windows(merged_peaks,1500)
dim(tfs_ar_1500)

