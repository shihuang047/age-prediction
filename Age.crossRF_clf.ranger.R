# ---
# title: "age prediction in the gut, oral and skin datasets"
# author: "ShiHuang"
# date: "12/26/2019"
# output: html_document
# ---
#-------------------------------
# install and load necessary libraries for data analyses
#-------------------------------

# install.packages('devtools') # if devtools not installed
# devtools::install_github('shihuang047/crossRanger')

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("biomformat")

p <- c("reshape2","ggplot2", "dplyr", "biomformat", "devtools", "crossRanger", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

#-------------------------------
# input args
#-------------------------------
setwd("/Users/huangshi/MyProjects/CMI-IBM/age-prediction/")
#-------------------------------
datafile<-"Input/gut_data/gut_4434.biom" # gut_data/gut_4434.biom | oral_data/oral_2550.biom | skin_data/skin_1975.biom
sample_metadata <- "Input/gut_data/gut_4434_map.txt" # gut_data/gut_4434_map.txt | oral_data/oral_2550_map.txt | skin_data/skin_1975_map.txt 
feature_metadata<-"Input/gut_data/gut_taxonomy.txt" # gut_data/gut_taxonomy.txt | oral_data/oral_taxonomy.txt | skin_data/skin_taxonomy.txt
prefix_name<-"gut_4434" # gut_4434 | oral_2550 | skin_1975
s_category<-c("cohort", "sex")  # c("cohort", "sex") | "qiita_host_sex" | c("body_site","qiita_host_sex") 
c_category<-"life_stage"  #"age" "qiita_host_age" "qiita_host_age"
outpath <- "./Output/gut_4434_by_cohort_sex_RF.clf_out/" # ./Output/gut_4434_by_cohort_sex_RF.reg_out/ ./Output/oral_2550_by_sex_RF.reg_out/ ./Output/skin_1975_by_site_sex_RF.reg_out/
dir.create(outpath)

#-------------------------------
# Biom table input
#-------------------------------
if(grepl("biom$", datafile)){
  biom <- read_biom(datafile)
  df <- data.frame(t(as.matrix(biom_data(biom))))
}else{
  df<-read.table(datafile, header=T, row.names=1, sep="\t", quote="", comment.char = "")
}

df<-df[order(rownames(df)), ]
 #df<-sweep(df, 1, rowSums(df), "/")
#-------------------------------
# Feature metadata input
#-------------------------------
fmetadata<-read.table(feature_metadata,header=T,sep="\t")
add_ann<-function(tab, fmetadata, tab_id_col=1, fmetadata_id_col=1){
  matched_idx<-which(fmetadata[, fmetadata_id_col] %in% tab[, tab_id_col])
  uniq_features_len<-length(unique(tab[, tab_id_col]))
  if(uniq_features_len>length(matched_idx)){
    warning("# of features has no matching IDs in the taxonomy file: ", uniq_features_len-length(matched_idx), "\n")
  }
  fmetadata_matched<-fmetadata[matched_idx,]
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
#-------------------------------
# Sample Metadata input
#-------------------------------
allmetadata<-read.table(sample_metadata,header=T,sep="\t",row.names=1, quote="", comment.char="")
if(length(allmetadata)==1){
  metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
  all_group<-colnames(metadata)<-colnames(allmetadata)
}else{
  metadata<-allmetadata[order(rownames(allmetadata)), ]
}
#-------------------------------
# Matching SampleID between biom data and metadata
#-------------------------------
## check if the order of rownames in the microbiome data and sample metadata are identical:
identical(rownames(df),rownames(metadata))
cat("The number of samples in biom table : ", nrow(df) ,"\n")
cat("The number of samples in metadata : ", nrow(metadata) ,"\n")

shared_idx<-intersect(rownames(df), rownames(metadata))
df_k<-df[shared_idx, ]
cat("The number of samples in biom (after filtering out samples with no metadata): ", nrow(df_k) ,"\n")
metadata_k<-metadata[shared_idx, ]
cat("The number of samples in metadata (after filtering out samples with no metadata): ", nrow(metadata) ,"\n")
identical(rownames(df_k),rownames(metadata_k))

#sink("./Input/oral_data/oral_2118_map.txt"); cat("#SampleID\t"); write.table(metadata_k, sep="\t", quote = F, row.names = T); sink()
#df<-df[,which(apply(df,2,var)!=0)]
#df_biom<-make_biom(as.matrix(t(df)))
#write_biom(df_biom, biom_file=file.path(paste("Input/gut_data/gut_", nrow(df),".biom", sep="")))
#-------------------------------
# To filter out samples with null values in both s and c categories
#-------------------------------
metadata_k<-metadata[which(!apply(metadata[, c(s_category, c_category)], 1,function(x) any(x==""))), ]
if(length(s_category)>1){
  metadata_k[, s_category]<-lapply(metadata_k[, s_category], factor)
}else{
  metadata_k[, s_category]<-factor(metadata_k[, s_category])
}
df_k<-df_k[rownames(metadata_k), ]
#-------------------------------
# fitler variables in microbiota data
#-------------------------------
cat("The number of samples : ", nrow(df_k) ,"\n")
cat("The number of variables : ", ncol(df_k) ,"\n")
#-------------------------------filtering taxa with zero variance
df_k<-df_k[,which(apply(df_k,2,var)!=0)]
cat("The number of fitlered variables (removed variables with zero variance) : ", ncol(df_k) ,"\n")
#-------------------------------filtering taxa with X% non-zero values
NonZero.p<-0.995
df_k<-df_k[,which(colSums(df_k==0)<NonZero.p*nrow(df_k))]
cat("The number of variables (removed variables containing over ", NonZero.p," zero) in training data: ", ncol(df_k) ,"\n")
#-------------------------------
# rf_clf using all datasets
#-------------------------------
x=df_k
y=factor(metadata_k[, c_category])

all_res_file<-paste(outpath, prefix_name, "_rf_clf_all_res.RData", sep="")
if(file.exists(all_res_file)){
  rf_all <- get(load(all_res_file))
}else{
  rf_all<-rf.cross.validation(x=x, y=y, ntree = 500, nfolds = 5, sparse=TRUE)
  save(rf_all, file=all_res_file)
}

#-------------------------------
# rf_reg with caret: tuning rf by 5-fold cv
#-------------------------------
library("caret")
#data<-data.frame(y, x)
# tgrid <- expand.grid(
#   .mtry = c(sqrt(ncol(x)), ncol(x)/3, ncol(x)), 
#   .splitrule = "gini",
#   .min.node.size = 5
# )
# tuned_rf <- train(y  ~ ., data = data,
#                      method = "ranger",
#                      trControl = trainControl(method="repeatedcv", number = 5, verboseIter = F),
#                      tuneGrid = tgrid,
#                      importance = 'impurity'
# )

#-------------------------------
# To filter out samples with null values in both s and c categories
#-------------------------------
metadata_k<-metadata[which(!apply(metadata[, c(s_category, c_category)], 1,function(x) any(x==""))), ]
if(length(s_category)>1){
  metadata_k[, s_category]<-lapply(metadata_k[, s_category], factor)
}else{
  metadata_k[, s_category]<-factor(metadata_k[, s_category])
}
df_k<-df_k[rownames(metadata_k), ]
#-------------------------------
# To creat a combined category if the length of s_category over 2
#-------------------------------
if(length(s_category)>=2){
  new_s_category<-paste0(s_category, collapse ="__")
  metadata[, new_s_category]<-do.call(paste, c(metadata[s_category], sep="_"))
  metadata_k[, new_s_category]<-do.call(paste, c(metadata_k[s_category], sep="_"))
  s_category=new_s_category
}
#-------------------------------
# rf_reg.by_datasets
#-------------------------------
## "rf_reg.by_datasets" runs standard random forests with oob estimation for regression of 
## c_category in each the sub-datasets splited by the s_category. 
## The output includes a summary of rf models in the sub datasets
## and all important statistics for each of features.
res_file<-paste(outpath, prefix_name, "_rf_clf.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_clf_res<-rf_clf.by_datasets(df_k, metadata_k, s_category, c_category,
                                 nfolds=3, verbose=FALSE, ntree=500)
  save(rf_clf_res, file=res_file)
}

rf_clf_res.summ<-plot_clf_res_list(rf_clf_res, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=outpath)

wilcox_res<-rf_clf_res.summ$feature_res
rf_models<-rf_clf_res$rf_model_list
# Add feature annotations using feature metadata
wilcox_res<-add_ann(wilcox_res, fmetadata)
q_cutoff=0.05
p.adj.method= "bonferroni"
## Non-specific features across datasets (at least present in two of datasets)
res_spcf<-id_non_spcf_markers(wilcox_res,positive_class="old", other_class="young", p.adj.method= "none", outdir=outpath)
summary(res_spcf)
wilcox_res_spcf<-res_spcf$feature_res_spcf
## keeps only OTUs which were significant in at least one datasets.
keep_sig_markers<-function(wilcox_res, q_cutoff){
  wilcox_res_sig<-wilcox_res[which(wilcox_res[, "feature"] %in% unique(as.character(subset(wilcox_res, Wilcoxon.test_p.adj < q_cutoff)[, "feature"]))), ]
  wilcox_res_sig
}

wilcox_res_spcf_sig<-keep_sig_markers(wilcox_res_spcf, q_cutoff=q_cutoff)
sink(paste(outpath,"BetweenGroupTest_out_all.xls",sep=""));write.table(wilcox_res_spcf,quote=FALSE,sep="\t", row.names = F);sink()
sink(paste(outpath,"BetweenGroupTest_out_sig_",p.adj.method,"_q_cutoff_",q_cutoff,".xls",sep=""));write.table(wilcox_res_spcf_sig,quote=FALSE,sep="\t", row.names = F);sink()


#-------------------------------
# rf_clf: subsampling at young ages 
# To address the reviewer's comment on why age prediction not accurate at old ages
#-------------------------------
ids_gt40<-which(metadata_k[, c_category]>40)
ids_st40<-which(metadata_k[, c_category]<40)
n_gt40<-length(ids_gt40)
n_st40<-length(ids_st40)

sub_ids_st40<-sample(ids_st40, n_gt40)
sub_x<-x[c(sub_ids_st40, ids_gt40),]; sub_x<-sub_x[,which(apply(sub_x,2,var)!=0)]; 
NonZero.p<-0.99
sub_x<-sub_x[,which(colSums(sub_x==0)<NonZero.p*nrow(sub_x))]
sub_y<-y[c(sub_ids_st40, ids_gt40)]
rf_sub<-rf.cross.validation(x=sub_x, y=sub_y, ntree = 500, nfolds = 5)

plot_obs_VS_pred(rf_sub$y, rf_sub$predicted, prefix="sub_train", target_field="age", span=1, outdir = outpath)
plot_perf_VS_rand(rf_sub$y, rf_sub$predicted, prefix="sub_train", target_field="age", n_features=ncol(sub_x), permutation = 1000, outdir = outpath)




