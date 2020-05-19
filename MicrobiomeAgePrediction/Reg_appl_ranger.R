#' ---
#' title: "age prediction/classification in the LS dataset"
#' author: "ShiHuang"
#' date: "8/6/2019"
#' output: html_document
#' ---
#'-------------------------------
#' install and load necessary libraries for data analyses
#'-------------------------------
p <- c("reshape2","ggplot2","pheatmap","combinat","plyr","ranger", "gridExtra",
       "vegan", "biomformat", "doMC", "cowplot", "pROC", "crossRanger")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#'-------------------------------
#' input args
#'-------------------------------
setwd("/Users/huangshi/MyProjects/CMI-IBM/Datasets/finrisk")
source("/Users/huangshi/MyProjects/CMI-IBM/R/data_trimming_util.R")
#'-------------------------------
train_datafile<-"gotu-tables/gotu.norm.biom"
train_sample_metadata<-"FR02.stool.txt"
train_feature_metadata<-NA
train_target_field<-"BL_AGE"
train_prefix<-"Finrisk"

test_datafile<-NA #gut_4575_rare_sp.csv
test_sample_metadata <- NA #"10283_20191126-092828.txt"
test_feature_metadata<-NA #"skin_taxonomy.txt"
test_prefix<-NA
test_target_field <- NA #

outpath <- "./Finrisk_gotu_RF.reg_out/"
h=40
AddTaxonomy=TRUE
dir.create(outpath)




#'-------------------------------
#' Train data: feature table input
#'-------------------------------
if(grepl("biom$", train_datafile)){
  train_biom <- read_biom(train_datafile)
  train_df <- data.frame(t(as.matrix(biom_data(train_biom))))
}else{
  train_df<-read.table(train_datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "")
}
train_df<-train_df[order(rownames(train_df)), ]
# df<-sweep(df, 1, rowSums(df), "/")
#'-------------------------------
#' Train data: sample Metadata input
#'-------------------------------
train_allmetadata<-read.table(train_sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
if(length(train_allmetadata)==1){
  train_metadata<-data.frame(train_allmetadata[order(rownames(train_allmetadata)),])
  train_all_group<-colnames(train_metadata)<-colnames(train_allmetadata)
}else{
  train_metadata<-train_allmetadata[order(rownames(train_allmetadata)), ]
  train_all_group<-colnames(train_metadata)
  train_all_group_f<-colnames(train_metadata)[sapply(train_metadata,class)=="factor"]
  train_all_group_n<-colnames(train_metadata)[sapply(train_metadata,class)!="factor"]
}

#'-------------------------------
#' Train data: filter out samples unmatched between biom and metadata
#'-------------------------------
train_data_list<-filter_samples_by_sample_ids_in_metadata(train_df, train_metadata)

#'-------------------------------
#' Train data: filter out samples with null values in train_target_field
#'-------------------------------
train_data_list<-filter_samples_by_NA_in_target_field_of_metadata(train_data_list$data, train_data_list$metadata, target_field = train_target_field)
#'-------------------------------
#' Train data: fitler features in feature table
#'-------------------------------
train_data_list$data<-filter_features_allzero(train_data_list$data)
cat("The number of kept features (removed features with zero variance) : ", ncol(train_data_list$data) ,"\n")
#-------------------------------filtering features with prevalence
prev=0.001
train_data_list$data<-filter_features_by_prev(train_data_list$data, prev= prev)
cat("The number of kept features (removed features containing less than ",  prev," prevalence) : ", ncol(train_data_list$data) ,"\n")

#'-------------------------------
#' trim the 16S feature length to 100nt 
#'-------------------------------
#test_data_list$data<-chop_seq_to_100nt(test_data_list$data)
#'-------------------------------
#' Train and test data: data input 
#'-------------------------------
train_x<-train_data_list$data
train_y<-train_data_list$metadata[, train_target_field]
train_y<-convert_y_to_numeric(train_y)


if(!is.na(test_datafile)){
  #'-------------------------------
  #' Test data: feature table input
  #'-------------------------------
  if(grepl("biom$", test_datafile)){
    test_biom <- read_biom(test_datafile)
    test_df <- data.frame(t(as.matrix(biom_data(test_biom))))
  }else{
    test_df<-read.table(datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "")
  }
  test_df<-test_df[order(rownames(test_df)),]
  #'df<-sweep(df, 1, rowSums(df), "/")
  #'-------------------------------
  #' Test data: feature metadata input
  #'-------------------------------
  if(!is.na(test_feature_metadata)){
    test_fmetadata<-read.table(test_feature_metadata,header=T,sep="\t")
  }
  #'-------------------------------
  #' Test data: sample Metadata input
  #'-------------------------------
  test_allmetadata<-read.table(test_sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
  if(length(test_allmetadata)==1){
    test_metadata<-data.frame(test_allmetadata[order(rownames(test_allmetadata)),])
    test_all_group<-colnames(test_metadata)<-colnames(test_allmetadata)
  }else{
    test_metadata<-test_allmetadata[order(rownames(test_allmetadata)), ]
    test_all_group<-colnames(test_metadata)
    test_all_group_f<-colnames(test_metadata)[sapply(test_metadata,class)=="factor"]
    test_all_group_n<-colnames(test_metadata)[sapply(test_metadata,class)!="factor"]
  }
  #'-------------------------------
  #' Test data: filter out samples unmatched between biom and metadata
  #'-------------------------------
  test_data_list<-filter_samples_by_sample_ids_in_metadata(test_df, test_metadata)
  #'-------------------------------
  #' Test data: filter out samples with null values in test_target_field
  #'-------------------------------
  test_data_list<-filter_samples_by_NA_in_target_field_of_metadata(test_data_list$data, test_data_list$metadata, target_field = test_target_field)
  #'-------------------------------
  #' Test data: filter out samples with too low read counts
  #'-------------------------------
  seq_dep_cutoff<-1000
  test_data_list<-filter_samples_by_seq_depth(test_data_list$data, test_data_list$metadata, cutoff = seq_dep_cutoff)
  #'-------------------------------
  #' Test data: fitler features by prevalence
  #'-------------------------------
  test_data_list$data<-filter_features_allzero(test_data_list$data)
  cat("The number of kept features (removed features with zero variance) : ", ncol(test_data_list$data) ,"\n")
  #'-------------------------------filtering features with prevalence
  prev=0.001
  test_data_list$data<-filter_features_by_prev(test_data_list$data, prev= prev)
  cat("The number of kept features (removed features containing less than ",  prev," prevalence) : ", ncol(test_data_list$data) ,"\n")
  #'-------------------------------
  #' Test data: rarefaction
  #'-------------------------------
  test_data_list$data<-rrarefy(test_data_list$data, seq_dep_cutoff)
  
  #'-------------------------------
  #' rf_reg in the test dataset itself
  #'-------------------------------
  x=test_data_list$data
  y=test_data_list$metadata[, test_target_field]
  # convert the y to a numberical variabley
  if(!is.numeric(y)){ y=as.numeric(as.character(y)) }
  res_file<-paste(outpath, test_prefix, "_rf_reg_res.RData", sep="")
  if(file.exists(res_file)){
    rf_reg_model<-get(load(res_file))
  }else{
    rf_reg_model<-rf.cross.validation(x, y, ntree = 5000, nfolds = 5)
    save(rf_reg_model, file=res_file)
  }
  plot_obs_VS_pred(rf_reg_model$y, rf_reg_model$predicted, target_field=paste(test_prefix, "age", sep=" "), span=1, outdir = outpath)
  plot_perf_VS_rand(rf_reg_model$y, rf_reg_model$predicted, target_field=paste(test_prefix, "age", sep=" "), permutation = 1000, n_features=NA, outdir = outpath)
  #plot_reg_feature_selection(x=x, y=y, rf_reg_model=rf_reg_model, metric ="MAE", outdir=outpath)
}

#'-------------------------------
#' Train and test data: keep features commonly shared btw train and test datasets
#'-------------------------------
if(!is.na(test_datafile)){
  print("yes")
  test_x<-test_data_list$data
  test_y<-test_data_list$metadata[, test_target_field]
  test_y<-convert_y_to_numeric(test_y)
  shared_data_list<-keep_shared_features(train_x, test_x)
  train_x_shared<-shared_data_list$train_x_shared
  test_x_shared<-shared_data_list$test_x_shared
}else{
  train_x_shared<-train_x
}

#'-------------------------------
#' Ranger modeling on train dataset
#'-------------------------------
train_res_file<-paste(outpath, train_prefix, "_shared_rf_reg_res.RData", sep="")
if(file.exists(train_res_file)){
  train_model<-get(load(train_res_file))
}else{
  train_model<-rf.out.of.bag(x=train_x_shared, y=train_y)
  save(train_model, file=train_res_file)
}
#'-------------------------------
#' Visualization of training model
#'-------------------------------
plot_obs_VS_pred(train_model$y, train_model$predicted, target_field="age", span=1, outdir = outpath)
plot_perf_VS_rand(train_model$y, train_model$predicted, target_field="age", permutation = 1000, n_features=NA, outdir = outpath)
plot_reg_feature_selection(x=train_x, y=train_y, rf_reg_model=train_model, metric ="MAE", outdir=outpath)

calc_rel_predicted(train_y, train_model$predicted, train_prefix="Finrisk", outdir = "./")



#'-------------------------------
#' Test data: rf modeling within test data using shared features
#'-------------------------------
test_shared_model<-rf.out.of.bag(x=test_x_shared, y=test_y)
plot_obs_VS_pred(test_shared_model$y, test_shared_model$predicted, prefix="test_shared", target_field="age", span=1, outdir = outpath)
plot_perf_VS_rand(test_shared_model$y, test_shared_model$predicted, prefix="test_shared", target_field="age", n_features = NA, permutation = 1000, outdir = outpath)
plot_residuals(test_shared_model$y, test_shared_model$predicted, prefix="test_shared", target_field="age", outdir = outpath)

#'-------------------------------
#' Test data: application of training model on test data using shared features
#'-------------------------------
test_pred<-predict(train_model$rf.model, test_x)
test_predicted<-test_pred$predictions
#'-------------------------------
#' Associate predicted y with metdata in the test data
#'-------------------------------
test_md_ma<-data.frame(test_data_list$metadata, test_predicted)
ggplot(test_md_ma, aes(x=get(test_target_field), y=test_predicted)) + geom_point()
#'-------------------------------
#' Test data: application of training model on test data using shared features
#'-------------------------------
plot_train_vs_test(train_y=train_y, predicted_train_y=train_model$predicted, 
                   test_y=test_y, predicted_test_y=test_pred$predictions, 
                   train_prefix=train_prefix, test_prefix=test_prefix, 
                   train_target_field="age", test_target_field="age", outdir=outpath)

calc_rel_predicted<-function(train_y, predicted_train_y, test_y=NULL, predicted_test_y=NULL,
                             train_prefix="train", test_prefix="test",
                             train_target_field="y", test_target_field=NULL, outdir=NULL){
  spl_train <- smooth.spline(train_y, predicted_train_y)
  train_relTrain <- residuals(spl_train); names(train_relTrain)<-names(train_y)
  train_fittedTrain <- fitted(spl_train); names(train_fittedTrain)<-names(train_y)
  relTrain_data<-train_relTrain_data <- data.frame(y=train_y, predicted_y=predicted_train_y, 
                                                   rel_predicted_y=train_relTrain, 
                                                   fitted_predicted_y=train_fittedTrain)
  sink(paste(outdir, train_prefix,".Relative_",train_target_field,".results.xls",sep=""));
  cat("\t");write.table(relTrain_data,quote=FALSE,sep="\t");sink()
  
  if(!is.null(test_y) & !is.null(predicted_test_y)){
    test_fittedTrain<-predict(spl_train, test_y)$y
    test_relTrain <- predicted_test_y - test_fittedTrain
    test_relTrain_data <- data.frame(y=test_y, predicted_y=predicted_test_y, 
                                     rel_predicted_y=test_relTrain,
                                     fitted_predicted_y=test_fittedTrain)
    relTrain_data<-rbind(train_relTrain_data, test_relTrain_data)
    DataSet <- factor(c(rep(train_prefix,length(train_relTrain)), rep(test_prefix,length(test_relTrain)) ))
    relTrain_data<-data.frame(relTrain_data, DataSet)
    sink(paste(outdir, train_prefix,"-",test_prefix,".Relative_",train_target_field,".results.xls",sep=""));
    cat("\t");
    write.table(relTrain_data,quote=FALSE,sep="\t");sink()
  }
  return(relTrain_data)
}
relTrain_data<-calc_rel_predicted(train_y=train_y, predicted_train_y=train_model$predicted, 
                                  test_y=test_y, predicted_test_y=test_predicted, 
                                  train_prefix=train_prefix, test_prefix=test_prefix, 
                                  train_target_field="age", test_target_field="age", outdir=outpath)

plot_residuals(test_y, test_predicted, prefix=test_prefix, target_field="age", outdir = outpath)

boxplot_rel_predicted_train_vs_test(relTrain_data, train_target_field="age", outdir=outpath)

plot_rel_predicted(relTrain_data, prefix = c(train_prefix, test_prefix), target_field = "age", outdir=outpath)

## plot train VS test data

p <- ggplot(relTrain_data, aes(x = y, y = predicted_y)) + 
  ylab(paste("Microbiome ", train_target_field, sep = "")) + 
  xlab(paste("Observed ", train_target_field, sep = "")) + 
  geom_point(aes(color = DataSet), alpha = 0.1) + 
  geom_smooth(aes(color = DataSet), method = "loess", span = 1)+
  facet_wrap(~DataSet) + theme(legend.position = "none") + 
  theme_minimal()
p
ggsave(filename = paste(outpath, train_prefix, "-", test_prefix, 
                        ".", test_target_field, ".train_test_ggplot.pdf", 
                        sep = ""), plot = p, height = 3, width = 6)
## boxplot of relative microbiome age
p_w <- formatC(stats::wilcox.test(rel_predicted_y ~ DataSet, data = relTrain_data)$p.value, digits = 4, format = "g")
p <- ggplot(relTrain_data, aes(x = DataSet, y = rel_predicted_y)) + 
  geom_violin(aes(color = DataSet)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4) + 
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.1) + 
  geom_hline(yintercept = 0) + theme_minimal() + 
  #ggtitle(paste("Wilcoxon Rank Sum Test:\n P=", p_w, sep = "")) + 
  xlab("Data sets") + 
  ylab(paste("Relative microbiome", train_target_field)) + 
  theme(legend.position = "none")
ggsave(filename = paste(outpath, train_prefix, "-", test_prefix, 
                        ".Relative_", train_target_field, ".boxplot.pdf", sep = ""), 
       width = 3, height = 4)
## age range in the test group
library("dplyr")
test_age_range<-range(relTrain_data[which(relTrain_data$DataSet==test_prefix),"y"])
relTrain_data_in_test_age_range<-relTrain_data %>% filter(y %in% (test_age_range[1]:test_age_range[2]))
plot_rel_predicted(relTrain_data_in_test_age_range, prefix = c(train_prefix, test_prefix), target_field = "age", outdir=outpath)



target_field="age"
prefix = c(train_prefix, test_prefix)
p <- ggplot(relTrain_data_in_test_age_range, aes(x = y, y = rel_predicted_y, color = DataSet)) + 
     ylab(paste("Relative microbiome ", target_field, sep = "")) + 
     xlab(paste("Observed ", target_field, sep = "")) + 
     geom_smooth(aes(group=DataSet),method = 'loess')+
     geom_point(alpha = 0.1) + 
     geom_hline(yintercept = 0) + theme_minimal()
p
prefix = paste(prefix, collapse = "-")
ggsave(filename = paste(outpath, prefix, ".", target_field, ".obs_vs_relative_pred.scatterplot.pdf", sep = ""), plot = p, device="pdf", 
       height = 4, width = 5)


## further see the daily changes of gut microbiome age
library(scales)
test_metadata$collection_time=as.Date(as.character(test_metadata$collection_timestamp), format = "%m/%d/%y")
test_data_list$metadata$collection_time=as.Date(as.character(test_data_list$metadata$collection_timestamp), format = "%m/%d/%y")
summary(relTrain_data_in_test_age_range)

## the selected metadata columns for further analysis.
selected_cols<-c("medical_intervention_time_automated", "resection_surgery_3_categories", "collection_time", "gluten.free", "vegan", "Green.Smoothies")
md_selected<-test_data_list$metadata[, selected_cols]
relTrain_data_in_test_age_range_md<-data.frame(subset(relTrain_data_in_test_age_range, DataSet=="Larry Smarr"), md_selected)
diet_summ<-data.frame(table(test_metadata[, c("vegan", "gluten.free", "Green.Smoothies")]))
diet_summ<-diet_summ[diet_summ[, "Freq"]>0,]
write.table(diet_summ, "diet_summ.xls", quote=F, sep="\t")
## plot: resection_surgery_3_categories
prefix = test_prefix
p <- ggplot(relTrain_data_in_test_age_range_md, 
            aes(x = collection_time, y = rel_predicted_y, color = resection_surgery_3_categories)) + 
  ylab(paste("Relative microbiome ", target_field, sep = "")) + 
  xlab("Date of sample collection") + 
  geom_smooth(method = 'loess')+
  geom_point(alpha = 0.1) + 
  #facet_wrap(~as.factor(y))+
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(0, 1), 
       legend.justification = c(0, 1), 
       legend.background = element_rect(colour = NA, fill = "white"))+
  theme_minimal()
  #theme_bw()
p

ggsave(filename = paste(outpath, prefix, ".", target_field, ".obs_vs_relative_pred.collection_time.scatterplot.pdf", sep = ""), plot = p, device="pdf", 
       height = 4, width = 10)
## plot: resection_surgery_3_categories + gluten.free
prefix = test_prefix
p <- ggplot(relTrain_data_in_test_age_range_md, 
            aes(x = collection_time, y = rel_predicted_y)) + 
  ylab(paste("Relative microbiome ", target_field, sep = "")) + 
  xlab("Date of sample collection") + 
  geom_smooth(method = 'loess')+
  geom_point(aes(color=gluten.free), alpha = 0.6) + 
  #facet_wrap(~as.factor(y))+
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1), 
        legend.background = element_rect(colour = NA, fill = "white"))+
  theme_minimal()
p
ggsave(filename = paste(outdir, prefix, ".", target_field, ".obs_vs_relative_pred.collection_time.gluten.free.scatterplot.pdf", sep = ""), plot = p, device="pdf", 
       height = 4, width = 10)

## plot: resection_surgery_3_categories + gluten.free
prefix = test_prefix
p <- ggplot(relTrain_data_in_test_age_range_md, 
            aes(x = collection_time, y = rel_predicted_y)) + 
  ylab(paste("Relative microbiome ", target_field, sep = "")) + 
  xlab("Date of sample collection") + 
  geom_smooth(method = 'loess')+
  #geom_line()+
  geom_point(aes(color=Green.Smoothies), alpha = 0.6) + 
  #facet_wrap(~as.factor(y))+
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1), 
        legend.background = element_rect(colour = NA, fill = "white"))+
  theme_minimal()
p
ggsave(filename = paste(outpath, prefix, ".", target_field, ".obs_vs_relative_pred.collection_time.Green.Smoothies.scatterplot.pdf", sep = ""), plot = p, device="pdf", 
       height = 4, width = 10)

## plot: resection_surgery_3_categories + gluten.free
prefix = test_prefix
p <- ggplot(relTrain_data_in_test_age_range_md, 
            aes(x = collection_time, y = rel_predicted_y)) + 
  ylab(paste("Relative microbiome ", target_field, sep = "")) + 
  xlab("Date of sample collection") + 
  geom_smooth(method = 'loess')+
  #geom_line()+
  geom_point(aes(color=vegan), alpha = 0.6) + 
  #facet_wrap(~as.factor(y))+
  geom_hline(yintercept = 0) + 
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1), 
        legend.background = element_rect(colour = NA, fill = "white"))+
  theme_minimal()
p
ggsave(filename = paste(outpath, prefix, ".", target_field, ".obs_vs_relative_pred.collection_time.vegan.scatterplot.pdf", sep = ""), plot = p, device="pdf", 
       height = 4, width = 10)


sessionInfo()
