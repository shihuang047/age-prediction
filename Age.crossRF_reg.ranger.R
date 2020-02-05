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
datafile<-"Input/oral_data/oral_2550.biom" # gut_data/gut_4434.biom | oral_data/oral_2550.biom | skin_data/skin_1975.biom
sample_metadata <- "Input/oral_data/oral_2550_map.txt" # gut_data/gut_4434_map.txt | oral_data/oral_2550_map.txt | skin_data/skin_1975_map.txt 
feature_metadata<-"Input/oral_data/oral_taxonomy.txt" # gut_data/gut_taxonomy.txt | oral_data/oral_taxonomy.txt | skin_data/skin_taxonomy.txt
prefix_name<-"oral_2550" # gut_4434 | oral_2550 | skin_1975
s_category<-"qiita_host_sex"  # c("cohort", "sex") | "qiita_host_sex" | c("body_site","qiita_host_sex") 
c_category<-"qiita_host_age"  #"age" "qiita_host_age" "qiita_host_age"
outpath <- "./Output/oral_2550_by_sex_RF.reg_out/" # ./Output/gut_4434_by_cohort_sex_RF.reg_out/ ./Output/oral_2550_by_sex_RF.reg_out/ ./Output/skin_1975_by_site_sex_RF.reg_out/
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
# rf_reg using all datasets
#-------------------------------
x=df_k
# convert the y to a numberical variable
if(!is.numeric(metadata[, c_category])){
  y=metadata_k[, c_category]=as.numeric(as.character(metadata_k[, c_category]))
}else{
  y=metadata_k[, c_category]
}

all_res_file<-paste(outpath, prefix_name, "_rf_reg_all_res.RData", sep="")
if(file.exists(all_res_file)){
  rf_all <- get(load(all_res_file))
}else{
  rf_all<-rf.cross.validation(x=x, y=y, ntree = 500, nfolds = 5, sparse=TRUE)
  save(rf_all, file=all_res_file)
}

cat(mean(rf_all$MAE), "+-", sd(rf_all$MAE), "\n")
cat(mean(rf_all$R_squared),"+-", sd(rf_all$R_squared),  "\n")
plot_obs_VS_pred(rf_all$y, rf_all$predicted, prefix="train", target_field="age", span=1, outdir = outpath)
# plot_perf_VS_rand(rf_all$y, rf_all$predicted, prefix="train", target_field="age", n_features=ncol(x), permutation = 1000, outdir = outpath)
#-------------------------------
# rf_reg with caret: tuning rf by 5-fold cv
#-------------------------------
library("caret")
data<-data.frame(y, x)
 tgrid <- expand.grid(
   .mtry = c(sqrt(ncol(x)), ncol(x)/3, ncol(x)), 
   .splitrule = "variance",
   .min.node.size = 5
 )
 tuned_rf <- train(y  ~ ., data = data,
                      method = "ranger",
                      trControl = trainControl(method="repeatedcv", number = 5, verboseIter = F),
                      tuneGrid = tgrid,
                      importance = 'impurity'
 )
 
## the default pars could be good enough

#rf.oob_all<-rf.out.of.bag(x=x, y=y, ntree = 500)
#plot_reg_feature_selection(x=x, y=rf.oob_all$y, rf.oob_all, outdir = outpath)
#-------------------------------
# groupKfold VS Kfold CV
# In this dataset, multiple samples may come from a same subject. 
# To test if the model can be generalizable across subject, we conducted the grouped K-fold CV, 
# where we kept samples of the same subject in either training or testing data.
#-------------------------------
data<-data.frame(y, x); dim(data)
# Kfold: run a random forest model in the regular K fold cv.
fit_control <- trainControl(## 5-fold CV
  method = "cv",
  number = 5)
set.seed(825)
rf_fit <- train(y ~ ., 
                data = data, 
                method = "ranger",
                trControl = fit_control)
rf_fit
# GroupKfold ensures samples from the same group (such as host subjects) do not present in both training and test dataset.
subjects<-metadata_k[, "host_subject_id"]
group_folds<-groupKFold(subjects, k=5)
group_fit_control <- trainControl(## use grouped CV folds
  index = group_folds,
  method = "cv",
  number=5)
set.seed(825)
rf_fit_group <- train(y ~ ., 
                data = data, 
                method = "ranger",
                trControl = group_fit_control)
rf_fit_group
save(rf_fit, file=paste(outpath, prefix_name, "-rf_reg_all_tuned_5fold_cv.rds", sep=""))
save(rf_fit_group, file=paste(outpath, prefix_name, "_rf_reg_all_tuned_group5fold_cv.rds", sep=""))
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
res_file<-paste(outpath, prefix_name, "_rf_reg.by_datasets_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_reg_res<-rf_reg.by_datasets(df_k, metadata_k, s_category, c_category,
                                 nfolds=3, verbose=FALSE, ntree=500)
  save(rf_reg_res, file=res_file)
}
## replace imp scores 0 with NA for features whose prevelance equal to 0 
for(i in 1:length(rf_reg_res$feature_imps_list)) {
  rf_reg_res$feature_imps_list[[i]][colSums(rf_reg_res$x_list[[i]])==0]<-NA
}
rf_reg_res$feature_imps_rank_list<-lapply(rf_reg_res$feature_imps_list, function(i) rank(-i, na.last = "keep"))
#rf_reg_res.summ<-plot.reg_res_list(rf_reg_res, outdir=outpath)
feature_res<-do.call(cbind, rf_reg_res$feature_imps_list)
feature_res_rank<-apply(feature_res, 2, function(x) rank(-x, na.last = "keep"))
rf_models<-rf_reg_res$rf_model_list
# Add feature annotations using feature metadata
feature_res<-data.frame(Feature.ID=rownames(feature_res), feature_res)
feature_res<-add_ann(feature_res, fmetadata)
feature_res_rank<-add_ann(data.frame(feature=rownames(feature_res_rank), feature_res_rank), fmetadata)
sink(paste(outpath,"feature_imps_all.xls",sep=""));write.table(feature_res,quote=FALSE,sep="\t", row.names = F);sink()
sink(paste(outpath,"feature_imps_rank_all.xls",sep=""));write.table(feature_res_rank,quote=FALSE,sep="\t", row.names = F);sink()

# Feature selection in all sub-datasets
# RF regression performance VS number of features used
# Prediction performances at increasing number of microbial species obtained by 
# retraining the random forest regressor on the top-ranking features identified 
# with a first random forest model training in a cross-validation setting
if(!file.exists(paste(outpath,"crossRF_feature_selection_summ.xls",sep=""))){
top_n_perf_list<-list()
for(n in 1:length(rf_reg_res$rf_model_list)){
  top_n_perf<-matrix(NA, ncol=5, nrow=11)
  max_n<-max(rf_reg_res$feature_imps_rank_list[[n]], na.rm = TRUE)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  colnames(top_n_perf)<-c("n_features", "MSE", "RMSE", "MAE", "MAE_perc")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  cat("Dataset: ", names(rf_reg_res$feature_imps_rank_list)[n], "\n")
  for(i in 1:length(n_features)){
    idx<-which(rf_reg_res$feature_imps_rank_list[[n]]<=n_features[i])
    x_n<-rf_reg_res$x_list[[n]][, idx]
    y_n<-rf_reg_res$y_list[[n]]
    top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500)
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-top_n_rf$MSE
    top_n_perf[i, 3]<-top_n_rf$RMSE
    top_n_perf[i, 4]<-top_n_rf$MAE
    top_n_perf[i, 5]<-top_n_rf$MAE_perc
  }
  top_n_perf[11, ]<-c(max_n, rf_reg_res$rf_MSE[n], rf_reg_res$rf_RMSE[n], rf_reg_res$rf_MAE[n], rf_reg_res$rf_MAE_perc[n])
  top_n_perf_list[[n]]<-top_n_perf
}
names(top_n_perf_list)<-names(rf_reg_res$rf_model_list)
top_n_perf_list<-lapply(1:length(top_n_perf_list), 
                        function(x) data.frame(Dataset=rep(names(top_n_perf_list)[x], nrow(top_n_perf_list[[x]])), top_n_perf_list[[x]]))
top_n_perf_comb<-do.call(rbind, top_n_perf_list)
top_n_perf_comb$n_features<-as.numeric(as.character(top_n_perf_comb$n_features))
top_n_perf_comb_m<-melt(top_n_perf_comb, id.vars = c("n_features", "Dataset"))
breaks<-top_n_perf_comb_m$n_features

p<-ggplot(subset(top_n_perf_comb_m, variable=="MAE"), aes(x=n_features, y=value)) + 
  xlab("# of features used")+
  ylab("MAE (yrs)")+
  scale_x_continuous(trans = "log",breaks=breaks)+
  geom_point(aes(color=Dataset)) + geom_line(aes(color=Dataset)) +#facet_wrap(~Dataset) +
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=18),
        strip.background = element_rect(colour = "white"), 
        panel.border = element_blank())
ggsave(filename=paste(outpath,"MAE__top_rankings.scatterplot.pdf",sep=""), plot=p, width=6, height=4)
sink(paste(outpath,"crossRF_feature_selection_summ.xls",sep=""));write.table(top_n_perf_comb,quote=FALSE,sep="\t", row.names = F);sink()
}

#-------------------------------
# rf_reg.cross_appl
#-------------------------------
# "rf_reg.cross_appl" runs standard random forests with oob estimation for regression of 
# c_category in each the sub-datasets splited by the s_category, 
# and apply the model to all the other datasets. 

crossRF_res<-rf_reg.cross_appl(rf_reg_res, rf_reg_res$x_list, rf_reg_res$y_list)
perf_summ<-crossRF_res$perf_summ
sink(paste(outpath,"crossRF_reg_perf_summ.xls",sep=""));write.table(perf_summ,quote=FALSE,sep="\t", row.names = F);sink()

#' The performance (MAE) of cross-applications  
#' The heatmap indicating MAE in the self-validation and cross-applications  
self_validation=as.factor(perf_summ$Train_data==perf_summ$Test_data)
library(viridis)
p_MAE<-ggplot(perf_summ, aes(x=as.factor(Test_data), y=as.factor(Train_data), z=MAE)) + 
  xlab("Test data")+ylab("Train data")+
  geom_tile(aes(fill = MAE, color = self_validation, width=0.9, height=0.9), size=1) + #
  scale_color_manual(values=c("white","grey80"))+
  geom_text(aes(label = round(MAE, 2)), color = "white") +
  scale_fill_viridis()+ 
  theme_bw() + theme_classic() +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = 90),
        axis.ticks = element_blank())
p_MAE
ggsave(filename=paste(outpath,"MAE_cross_appl_matrix_",c_category, "_among_", s_category,".heatmap.pdf",sep=""),plot=p_MAE, width=5, height=4)

#' The scatter plot matrix showing predicted and Reported values in the self-validation and cross-applications 
predicted_summ<-dplyr::bind_rows(crossRF_res$predicted, .id = "Train_data__VS__test_data")
tmp<-data.frame(do.call(rbind, strsplit(predicted_summ$Train_data__VS__test_data, "__VS__")))
colnames(tmp)<-c("Train_data", "Test_data")
self_validation=as.factor(tmp$Train_data==tmp$Test_data)
predicted_summ<-data.frame(tmp, self_validation, predicted_summ)

#l<-levels(data$sex); l_sorted<-sort(levels(data$sex))
Mycolor <- c("#0072B2", "#D55E00") 
#if(identical(order(l), order(l_sorted))){Mycolor=Mycolor }else{Mycolor=rev(Mycolor)}
target_variable="age"
p_scatter<-ggplot(predicted_summ, aes(x=test_y, y=pred_y))+
         ylab(paste("Predicted ",target_variable,sep=""))+
         xlab(paste("Reported ",target_variable,sep=""))+
         geom_point(aes(color=self_validation), alpha=0.1)+
         geom_smooth(aes(color=self_validation), method="loess",span=1)+
         scale_color_manual(values = Mycolor)+
         facet_grid(Train_data~Test_data)+
        theme_bw()+
        theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank())+
        theme(legend.position="none")
p_scatter
ggsave(filename=paste(outpath,"Scatterplot_cross_appl_matrix_",c_category, "_among_", s_category,".pdf",sep=""),plot=p_scatter, width=8, height=8)

if(grepl("skin", prefix_name)){
# Test if skin microbiota age prediction model can be applied across forehead and hand microbiota
s_category="body_site"
res_file<-paste(outpath, prefix_name, "_rf_reg.by_", s_category ,"_res.RData", sep="")
if(file.exists(res_file)){
  load(res_file)
}else{
  rf_reg_res<-rf_reg.by_datasets(df_k, metadata_k, s_category, c_category,
                                 nfolds=3, verbose=FALSE, ntree=500)
  save(rf_reg_res, file=res_file)
}

crossRF_res<-rf_reg.cross_appl(rf_reg_res, rf_reg_res$x_list, rf_reg_res$y_list)
perf_summ<-crossRF_res$perf_summ
sink(paste(outpath, s_category, "crossRF_reg_perf_summ.xls",sep=""));write.table(perf_summ,quote=FALSE,sep="\t", row.names = F);sink()


#' The scatter plot matrix showing predicted and Reported values in the self-validation and cross-applications 
predicted_summ<-dplyr::bind_rows(crossRF_res$predicted, .id = "Train_data__VS__test_data")
tmp<-data.frame(do.call(rbind, strsplit(predicted_summ$Train_data__VS__test_data, "__VS__")))
colnames(tmp)<-c("Train_data", "Test_data")
self_validation=as.factor(tmp$Train_data==tmp$Test_data)
predicted_summ<-data.frame(tmp, self_validation, predicted_summ)
#l<-levels(data$sex); l_sorted<-sort(levels(data$sex))
Mycolor <- c("#0072B2", "#D55E00") 
#if(identical(order(l), order(l_sorted))){Mycolor=Mycolor }else{Mycolor=rev(Mycolor)}
target_variable="age"
p_scatter<-ggplot(predicted_summ, aes(x=test_y, y=pred_y))+
  ylab(paste("Predicted ",target_variable,sep=""))+
  xlab(paste("Reported ",target_variable,sep=""))+
  geom_point(aes(color=Test_data), alpha=0.1)+
  geom_smooth(aes(color=Test_data), method="loess",span=1)+
  scale_color_manual(values = Mycolor)+
  facet_grid(Train_data~Test_data)+
  theme_bw()+
  theme(axis.line = element_line(color="black"),
        strip.background = element_rect(colour = "white"),
        panel.border = element_blank())+
  theme(legend.position="none")
p_scatter
ggsave(filename=paste(outpath,"Scatterplot_cross_appl_matrix_",c_category, "_among_", s_category,".pdf",sep=""),plot=p_scatter, width=4, height=4)


# The histograms showing the prediction accuracy of age regression models 
# dependent on skin body sites and their cross-applications compared to random permutations.
  target_variable="age"
  rand_MAE <- function(y, permutation=1000) {
    set.seed(123)
    rand_y_mat <- replicate(permutation, sample(y, replace = FALSE))
    apply(rand_y_mat, 2, function(x) MAE <- mean(sqrt((y - x)^2)) )
  }
  rand_MAE_list<-lapply(crossRF_res$predicted, function(x) data.frame(rand_MAE=rand_MAE(x$test_y), MAE=mean(sqrt((x$test_y - x$pred_y)^2))))
  rand_MAE_df<-do.call(rbind, rand_MAE_list)
  rand_MAE_df$Dataset<-rep(names(crossRF_res$predicted), each=1000)
  rand_MAE_df$Train_data<-rep(perf_summ$Train_data, each=1000)
  rand_MAE_df$Test_data<-rep(perf_summ$Test_data, each=1000)
  
  p_hist<-ggplot(rand_MAE_df, aes(x = rand_MAE)) + 
    geom_histogram(alpha = 0.5) + 
    xlab("MAE") + 
    geom_vline(data = rand_MAE_df, aes(color=Test_data, xintercept = MAE)) + 
    scale_color_manual(values = Mycolor)+
    facet_grid(Train_data~Test_data)+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())+
    theme(legend.position="none")
  p_hist
  ggsave(filename=paste(outpath,"Histogram_cross_appl_rand_MAE_matrix_",c_category, "_among_", s_category,".pdf",sep=""),plot=p_hist, width=4, height=4)
}

#-------------------------------
# rf_reg: subsampling at young ages 
# To address the reviewer's comment on why age prediction not accurate at old ages
#-------------------------------
ids_gt40<-which(metadata_k[, c_category]>40)
ids_st40<-which(metadata_k[, c_category]<=40)
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




