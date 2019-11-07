#--------------------------------------------------
p <- c("reshape2","randomForest", "optparse", "ade4", "doMC",
       "ggplot2", "RColorBrewer", "vegan", "xgboost", "caret")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))

normalize_NA_in_metadata<-function(md){
  apply(md, 1, function(x) {
   idx <- which(x=="not provided" | x=="Not provided" | x=="Not Provided"
          | x=="not applicable" | x=="Not applicable"
          | x=="NA" | x=="na" | x=="Na"
          | x=="none" | x=="None" | x=="NONE")
   x[idx]<-NA
  })
  md
}
discard_uninfo_columns_in_metadata<-function(md){
  noninfo_idx<-which(apply(md, 2, function(x) length(unique(x))==1))
  md<-md[-noninfo_idx]
  md
}


trim_metadata<-function(md){
  if(length(md)==1){
    md<-data.frame(md[order(rownames(md)),])
    test_all_group<-colnames(md)<-colnames(md)
  }else{
    md<-md[order(rownames(md)), ]
    md<-normalize_NA_in_metadata(md)
    md<-discard_uninfo_columns_in_metadata(md)
    #all_group<-colnames(md)
    #all_group_f<-colnames(md)[sapply(md,class)=="factor"]
    #all_group_n<-colnames(md)[sapply(md,class)!="factor"]
  }
  md
}

filter_features_allzero<-function(data, samples=TRUE, features=TRUE){
  if(samples & features){
    result<-data[which(apply(data, 1, sum)!=0), ]
    result<-data[, which(apply(result, 2, sum)!=0)]
  }else if(samples & !features){
    result<-data[which(apply(data, 1, sum)!=0), ]
  }else if(!samples & features){
    result<-data[, which(apply(data, 2, sum)!=0)]
  }else{
    stop("Nothing has been done!")
  }
  result
}

filter_features_by_prev <- function(data, prev=0.001){
  data<-data[, which(colSums(data!=0) > prev * nrow(data))]
  data
}


filter_samples_by_NA_in_y <- function(data, y){
  y_k<-y[which(!is.na(y))]
  data_k<-data[which(!is.na(y)) ,]
  result<-list()
  result$data_k<-data_k
  result$y_k<-y_k
  result
}

filter_samples_by_NA_in_target_field_of_metadata <- function(data, metadata, target_field){
  if(!identical(rownames(data), rownames(metadata)))stop("The sample IDs should be idenical in feature table and metadata!")
  idx<-which(metadata[, target_field]!='not applicable' & 
                  metadata[, target_field]!="not provided" & 
                  !is.na(metadata[, target_field]) )
  metadata_k<-metadata[idx, ]
  data_k<-data[idx, ]
  cat("The number of samples kept (after filtering out samples with NA values in ",target_field,"): ", nrow(metadata_k) ,"\n")
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  result
}


filter_samples_by_sample_ids_in_metadata <- function(data, metadata){
  shared_idx<-intersect(rownames(data), rownames(metadata))
  data_matched<-data[shared_idx, ]
  cat("The number of samples in feature table (after filtering out samples with no metadata): ", nrow(data_matched) ,"\n")
  metadata_matched<-metadata[shared_idx, ]
  cat("The number of samples metadata (after filtering out samples with no metadata): ", nrow(metadata_matched) ,"\n")
  cat("The sample IDs are idenical in feature table and metadata: ", identical(rownames(data_matched), rownames(metadata_matched)), "\n")
  result<-list()
  result$data<-data_matched
  result$metadata<-metadata_matched
  return(result)
}

filter_samples_by_seq_depth<-function(data, metadata, cutoff=1000){
  seq_dep_idx<-which(rowSums(data) > cutoff)
  cat("The number of kept samples with more than ",cutoff," reads: ", length(seq_dep_idx), "\n")
  if(length(seq_dep_idx)>0){
    metadata_k<-metadata[seq_dep_idx, ]
    data_k<-data[seq_dep_idx, ]
  }else{
    stop("The read count of all samples is less than sequencing depth threshold!")
  }
  
  #metadata_k$Seq_depth<-rowSums(df_k)
  #p<-ggplot(metadata_k, aes(Seq_depth)) + geom_histogram() + xlim(1000, 70000) + theme_bw()
  #p
  result<-list()
  result$data<-data_k
  result$metadata<-metadata_k
  return(result)
}

convert_y_to_numeric<-function(y, reg=TRUE){
  if(reg & !is.numeric(y)){
    y=as.numeric(as.character(y))
  }else if(!reg){
    y=factor(y)
  }else{
    y=y
  }
  y
}

keep_shared_features<-function(train_x, test_x){
  common_idx<-intersect(colnames(train_x), colnames(test_x))
  train_x_shared<-train_x[, common_idx]
  test_x_shared<-test_x[, common_idx]
  cat("Number of features kept:", length(common_idx), "\n")
  cat("The proportion of commonly shared features in train and test dataset respectively: \n")
  cat("Train data: ", length(common_idx)/ncol(train_x), "\n")
  cat("Test data: ", length(common_idx)/ncol(test_x), "\n")
  result<-list()
  result$train_x_shared<-train_x_shared
  result$test_x_shared<-test_x_shared
  result
}


#' feature metadata
add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  fmetadata_matched<-fmetadata[which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col]),]
  out<-merge(tab, fmetadata_matched, by.x=tab_id_col, by.y=fmetadata_id_col)
  out
}
rbind.na<-function(l){
  max_len<-max(unlist(lapply(l, length)))
  c_l<-lapply(l, function(x) {c(x, rep(NA, max_len - length(x)))})
  do.call(rbind, c_l)
}
expand_Taxon<-function(df, Taxon){
  taxa_df <- rbind.na(strsplit(as.character(df[, Taxon]), '; '))
  colnames(taxa_df) <- c("kingdom","phylum","class","order","family","genus","species") #"kingdom", 
  data.frame(df, taxa_df)
}
