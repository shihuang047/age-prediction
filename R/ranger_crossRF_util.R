#' Loading libraries

p <- c("ade4","caret","reshape2","MASS","ranger","plyr", "foreach", "doMC", "ALEPlot", "vegan", "compositions",
       "microbenchmark", "lattice","gplots","ggplot2","squash","RColorBrewer","pheatmap","pROC","ROCR","caTools")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))


#' customised colors for Enr's 3 factors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#' comb function for parallelization using foreach
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#' Perform pairwise rf classfication for a data matrix between all pairs of group levels
#' @param x a data matrix or data.frame
#' @param i the order of one level
#' @param j the order of the other level
#' @param f factor
#' @return result list of the random forest classification
#'
#' @seealso ranger
#' @examples
#' 
#' df <- t(rmultinom(16,160,c(.001,.6,.2,.3,.299))) + 0.65
#' f<-factor(c(rep("A", 4), rep("B", 4), rep("C", 4), rep("D", 4)))
#' f0<-factor(c(rep("A", 8), rep("B", 8)))
#' rf_compare_levels(df=df, f=f, i=1, j=2)
#' rf_clf.pairwise(df, f0)
#' 

rf_clf.pairwise <- function (df, f, nfolds=3, ntree=5000, verbose=FALSE) {
  rf_compare_levels <- function(df, f, i=1, j=2, nfolds=3, ntree=500, verbose=FALSE) {
    df_ij <- df[which(as.integer(f) == i | as.integer(f) == j), ]
    f_ij <- factor(f[which(as.integer(f) == i | as.integer(f) == j)])
    if(nfolds==3){
      oob<-rf.out.of.bag(df_ij, f_ij, verbose=verbose, ntree=ntree)
    }else{
      oob<-rf.cross.validation(df_ij, f_ij, nfolds=nfolds, verbose=verbose, ntree=ntree)
    }
    positive_class=levels(f)[i]
    cat("\nTraining dataset: ", levels(f)[i], "-", levels(f)[j] ,"\n\n")
    conf<-caret::confusionMatrix(data=oob$predicted, f_ij, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the cross-validation: ", acc ,"\n") 
    rf_AUC<-get.auroc(oob$probabilities, f_ij, positive_class)
    cat("AUC in the cross-validation: ", rf_AUC ,"\n") # wired value of "1" kept showing
    c("AUC"=rf_AUC, acc, kappa_oob, conf$byClass)
  }
  if(nlevels(f)==2){
    out_summ<-rf_compare_levels(df, f, i=1, j=2, nfolds, ntree, verbose)
  }else{
    level_names<-levels(f)
    ix <- setNames(seq_along(level_names), level_names)
    out_list<-outer(ix[-1L], ix[-length(ix)], 
                    function(ivec, jvec) sapply(seq_along(ivec), 
                                                function(k) {
                                                  i <- ivec[k]
                                                  j <- jvec[k]
                                                  if (i > j) {
                                                    rf_compare_levels(df, f, i, j, nfolds, ntree, verbose)
                                                  }else{NA}
                                                }))
    dataset_name<-outer(dimnames(out_list)[[1]],dimnames(out_list)[[2]], paste, sep="__")
    out_rownames<-dataset_name[lower.tri(dataset_name, diag = T)]
    out_summ<-do.call(rbind, out_list[lower.tri(out_list, diag = T)]); rownames(out_summ)<-out_rownames
  }
  out_summ
}


#' Runs standard random forests with oob estimation for classification of 
#' c_category in each the sub-datasets splited by the s_category. 
#' The output includes a summary of rf models in the sub datasets
#' and all important statistics for each of features.
rf_clf.by_datasets.summ<-function(df, metadata, s_category, c_category, positive_class=NA, nfolds=3, verbose=FALSE, ntree=5000, 
                                  p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir="./"){
  res_list<-rf_clf.by_datasets(df, metadata, s_category, c_category, positive_class, nfolds, verbose, ntree)
  stopifnot(all(names(res_list) %in% c("datasets","rf_model_list","sample_size","rf_AUC","feature_imps_list")==TRUE))
  plot_res_list<-plot.res_list(res_list, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=outdir)
  result<-list()
  result$rf_models<-res_list$rf_model_list
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  return(result)
}


#' Runs standard random forests with oob estimation for classification of 
#' c_category in each the sub-datasets splited by the s_category, 
#' and apply the model to all the other datasets. The output includes
#' accuracy, auc and Kappa statistics.
#'
#' @param df Training data: a data.frame.
#' @param metadata Sample metadata with at least two columns.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'factor' used as sample label for rf classification in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. 
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))), 
#'                      f_c=factor(c(rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8))),
#'                      age=c(rep("C", 15), rep("H", 15), rep("C", 15), rep("H", 15)))
#'                      
#' 
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c', positive_class="C")
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c', positive_class="C", rf_imp_pvalues=TRUE)
#'
#'
rf_clf.by_datasets<-function(df, metadata, s_category, c_category, positive_class=NA, rf_imp_pvalues=FALSE, clr_transform=TRUE, nfolds=3, verbose=FALSE, ntree=500,
                             p.adj.method = "BH", q_cutoff=0.05){
  y_list<-split(metadata[, c_category], metadata[, s_category])
  x_list<-split(df, metadata[, s_category])
  datasets<-levels(factor(metadata[, s_category]))
  L<-length(y_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  # 1. sample size of all datasets
  sample_size<-as.numeric(table(metadata[, s_category]))
  nCores <- detectCores()
  registerDoMC(nCores)
  oper<-foreach(i=1:L, .combine='comb', .multicombine=TRUE, .init=list(list(), list(), list())) %dopar% {
    x<-x_list[[i]]
    y<-factor(y_list[[i]])
    # 2. AUC of random forest model
    if(nfolds==3){
      oob<-rf.out.of.bag(x, y, verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
      rf_imps=oob$importances
    }else{
      oob<-rf.cross.validation(x, y, nfolds=nfolds, verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
      rf_imps=rowMeans(oob$importances)
    }
    rf_AUC<-get.auroc(oob$probabilities, y, positive_class)
    # 3. # of significantly differential abundant features between health and disease
    out<-BetweenGroup.test(x, y, clr_transform=clr_transform, positive_class=positive_class, p.adj.method = p.adj.method, q_cutoff=q_cutoff)
    feature_imps<-data.frame(feature=rownames(out), dataset=rep(datasets[i], ncol(x)), 
                             rf_imps=rf_imps, out)
    list(oob=oob, rf_AUC=rf_AUC, feature_imps=feature_imps)
  }
  
  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-oper[[1]]
  result$rf_AUC<-unlist(oper[[2]])
  result$feature_imps_list<-oper[[3]]
  
  return(result)
}

#' Runs standard random forests with oob estimation for regression of 
#' c_category in each the sub-datasets splited by the s_category, 
#' and apply the model to all the other datasets. The output includes
#' accuracy, auc and Kappa statistics.
#'
#' @param df Training data: a data.frame.
#' @param metadata Sample metadata with at least two columns.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'factor' used as sample label for rf classification in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. 
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_s1=factor(c(rep(TRUE, 30), rep(FALSE, 30))),
#'                      f_c=factor(c(rep("C", 15), rep("H", 15), rep("D", 15), rep("P", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' 
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_s', c_category='age')
#' reg_res
rf_reg.by_datasets<-function(df, metadata, s_category, c_category, nfolds=3, rf_imp_pvalues=FALSE, verbose=FALSE, ntree=500){
  #as.numeric.factor <- function(y) {as.numeric(levels(y))[y]}
  y<-metadata[, c_category]
  if(is.factor(metadata[, c_category])) y<-as.numeric(as.character(y))
  y_list<-split(y, metadata[, s_category])
  x_list<-split(df, metadata[, s_category])
  # data check
  try(if(!all(unlist(lapply(y_list, length)) == unlist(lapply(x_list, nrow)))) 
    stop("# of samples in x and length of y should match to each other within all datasets!") )
  datasets<-levels(factor(metadata[, s_category]))
  L<-length(y_list)
  # 1. sample size of all datasets
  sample_size<-as.numeric(table(metadata[, s_category]))
  nCores <- detectCores()
  registerDoMC(nCores)
  oper<-foreach(i=1:L, .combine='comb', .multicombine=TRUE, 
                .init=list(list(), list(), list(), list(), 
                           list(), list(), list(), list(), list())) %dopar% {
                             # 2. AUC of random forest model
                             if(nfolds==3){
                               oob<-rf.out.of.bag(x_list[[i]], y_list[[i]], verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
                               rf_imps=oob$importances
                               rf_model<-oob$rf.model
                             }else{
                               oob<-rf.cross.validation(x_list[[i]], y_list[[i]], nfolds=nfolds, verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
                               rf_model<-oob$models
                               rf_imps=rowMeans(oob$importances)
                             }
                             MSE=oob$MSE
                             RMSE=oob$RMSE
                             MAE=oob$MAE
                             MAE_perc=oob$MAE_perc
                             R_squared=oob$R_squared
                             Adj_R_squared=oob$Adj_R_squared
                             predicted=oob$predicted
                             list(rf_model=rf_model, predicted=predicted, MSE=MSE, RMSE=RMSE,
                                  MAE=MAE, MAE_perc=MAE_perc, R_squared=R_squared, Adj_R_squared=Adj_R_squared, feature_imps=rf_imps)
                           }
  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-oper[[1]]; names(result$rf_model_list)<-result$datasets
  result$rf_predicted<-oper[[2]]; names(result$rf_predicted)<-result$datasets
  result$feature_imps_list<-oper[[9]]; names(result$feature_imps_list)<-result$datasets
  result$rf_MSE<-unlist(oper[[3]])
  result$rf_RMSE<-unlist(oper[[4]])
  result$rf_MAE<-unlist(oper[[5]])
  result$rf_MAE_perc<-unlist(oper[[6]])
  result$rf_R_squared<-unlist(oper[[7]])
  result$rf_Adj_R_squared<-unlist(oper[[8]])
  return(result)
}

#' The correlation of importance scores between datasets 
#' @param reg_res a list, the output of the function 'rf_reg.by_datasets'
#' @param ranked if transform importance scores into rank for correlation analysis
#' @param plot if plot the correlation matrix
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_c=factor(c(rep("C", 15), rep("H", 15), rep("D", 15), rep("P", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' 
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_c', c_category='age')
#' corr_datasets_by_imps(reg_res$feature_imps_list)
#' 

corr_datasets_by_imps<-function(feature_imps_list, ranked=TRUE, plot=FALSE){
  if(ranked){
    ranked_list<-lapply(feature_imps_list, function(x) rank(-x))
    corr_mat<-cor(do.call(cbind, ranked_list))
  }else{
    corr_mat<-cor(do.call(cbind, feature_imps_list))
  }
  if(plot){
    require(corrplot)
    corrplot(corr =corr_mat, order="AOE", type="upper") # tl.pos = "d"
    corrplot(corr = corr_mat, add=TRUE, type="lower", method="number",
             order="AOE",diag=TRUE, tl.pos="n", cl.pos="n")
  }
  corr_mat
}


plot.clf_res_list<-function(res_list, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir="./"){
  datasets<-res_list$datasets
  sample_size<-res_list$sample_size
  rf_AUC<-res_list$rf_AUC
  feature_imps_list<-res_list$feature_imps_list
  # Statistics summary of biomarkers discovery in multiple datasets
  num_sig_p.adj<-sapply(feature_imps_list, function(x) sum(x[,"Wilcoxon.test_p.adj"] < q_cutoff))
  num_sig_p<-sapply(feature_imps_list, function(x) sum(x[,"Wilcoxon.test_p"]< p_cutoff))
  num_enriched<-sapply(feature_imps_list, function(x) length(grep("enriched", x[,"Enr"])))
  num_depleted<-sapply(feature_imps_list, function(x) length(grep("depleted", x[,"Enr"])))
  # ggplot
  summ<-data.frame(datasets=datasets, sample_size=sample_size, AUC=rf_AUC, num_sig_p, num_sig_p.adj, num_enriched, num_depleted)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  #summ_m<-melt(summ)
  p_a<- ggplot(summ, aes(x=Data_sets, y=Sample_size)) + xlab("Data sets") + ylab("Sample size")+ 
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=Data_sets, y=AUC)) + xlab("") + ylab("AUROC")+ 
    geom_point(shape="diamond", size=4)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    geom_hline(yintercept=0.5, linetype="dashed")+
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    scale_y_continuous(limits = c(0.5,1)) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=Data_sets, y=num_sig_p.adj)) + xlab("") + ylab("# of sig. features")+ 
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  summ_Enr<-melt(summ[, c("Data_sets","num_enriched","num_depleted")])
  summ_Enr[summ_Enr$variable=="num_depleted",]$value<--summ_Enr[summ_Enr$variable=="num_depleted",]$value
  p_d<- ggplot(summ_Enr, aes(x=Data_sets, y=value, colour=variable)) + xlab("") + ylab("Enrichment")+ 
    ylim(-max(summ_Enr$value), max(summ_Enr$value))+
    geom_hline(yintercept=0)+
    geom_point(show.legend=F, size=2)+ 
    geom_bar(show.legend=F, stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  require('gridExtra')
  p1<-arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(4, 2, 2, 2))
  ggsave(filename=paste(outdir,"Datasets_AUC.ggplot.pdf",sep=""),p1, width=9, height=3+nrow(summ)*0.2)
  # boxplot indicating p and p.adj values of sig. features
  feature_res<-plyr::ldply(feature_imps_list)
  feature_res_m<-melt(feature_res[, c("feature","dataset", "Enr", "AUC","rf_imps", "Wilcoxon.test_p", "Wilcoxon.test_p.adj")])
  # customised colors for Enr's 3 factors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  my3cols<-c("grey60", rev(gg_color_hue(2)))
  p2<-ggplot(feature_res_m, aes(x=dataset, y=value)) + 
    geom_boxplot(outlier.shape = NA)+
    facet_grid(.~variable, scales="free") + scale_color_manual(values = my3cols) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    geom_jitter(aes(color=Enr), position=position_jitter() ,size=1, alpha=0.4)+ #jitter
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"), 
          panel.border = element_blank()) 
  ggsave(filename=paste(outdir,"Datasets_feature_res.ggplot.pdf",sep=""),plot=p2, width=10, height=4)
  # library(ggpubr)
  # ggarrange(plotlist = list(p1, p2))
  result<-list()
  result$summ<-summ
  result$summ_plot<-p1
  result$feature_res<-feature_res
  result$feature_res_plot<-p2
  return(result)
}


plot.reg_res_list<-function(res_list, outdir="./"){
  datasets<-res_list$datasets
  sample_size<-res_list$sample_size
  rf_MSE<-res_list$rf_MSE
  rf_RMSE<-res_list$rf_RMSE
  rf_MAE<-res_list$rf_MAE
  rf_MAE_perc<-res_list$rf_MAE_perc
  rf_R_squared<-res_list$rf_R_squared
  rf_Adj_R_squared<-res_list$rf_Adj_R_squared
  feature_imps_list<-res_list$feature_imps_list
  # prevelance of features in subdatasets
  prev_df<-do.call(cbind, lapply(res_list$x_list, function(x) apply(x, 2, function(a) sum(a==0)/length(a))))
  # ggplot
  summ<-data.frame(datasets=datasets, sample_size=sample_size, 
                   MSE=rf_MSE, RMSE=rf_RMSE, MAE=rf_MAE, MAE_perc=rf_MAE_perc, 
                   R_squared=rf_R_squared, Adj_R_squared=rf_Adj_R_squared)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  #summ_m<-melt(summ)
  p_a<- ggplot(summ, aes(x=Data_sets, y=Sample_size)) + xlab("Data sets") + ylab("Sample size")+ 
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=Data_sets, y=RMSE)) + xlab("") + ylab("RMSE")+ 
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=Data_sets, y=MAE)) + xlab("") + ylab("MAE")+ 
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_d<- ggplot(summ, aes(x=Data_sets, y=R_squared)) + xlab("") + ylab("R_squared")+ 
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  require('gridExtra')
  summ_plot<-arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(3, 2, 2, 2))
  ggsave(filename=paste(outdir,"Datasets_perfs.ggplot.pdf",sep=""), summ_plot, width=9, height=3+nrow(summ)*0.2)
  # boxplot indicating p and p.adj values of sig. features
  #feature_res<-plyr::ldply(feature_imps_list)
  feature_res<-do.call(cbind, feature_imps_list)
  result<-list()
  result$summ<-summ
  result$summ_plot<-summ_plot
  result$feature_res<-feature_res
  return(result)
}

#' Based on pre-computed rf models classifying 'c_category' in each the sub-datasets splited by the 's_category',
#' perform cross-datasets application of the rf models. The inputs are precalculated 
#' rf models, and the outputs include accuracy, auc and Kappa statistics.
#' 
#' @param rf_model_list A list of rf.models generated from the function rf.out.of.bag.
#' @param x_list A list of training data: data.frame.
#' @param y_list A list of y.
#' @param positive_class A string indicates one common class in each of elements in the y_list.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      f_s0=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_c=factor(c(rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8))),
#'                      age=c(1:30, 2:31)
#'                      )
#' res_list<-rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c', positive_class="C")
#' rf_model_list<-res_list$rf_model_list
#' rf_clf.cross_appl(rf_model_list, res_list$x_list, res_list$y_list, positive_class="C")
#' 
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f=metadata[, 'f_s'], comp_group, verbose=FALSE, ntree=500, p.adj.method = "bonferroni", q_cutoff=0.05)                
#' comps_res 
#' rf_clf.cross_appl(comps_res$rf_model_list, x_list=comps_res$x_list, y_list=comps_res$y_list, positive_class=comp_group)
#'  
rf_clf.cross_appl<-function(rf_model_list, x_list, y_list, positive_class=NA){
  L<-length(rf_model_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  try(if(!identical(L, length(x_list), length(y_list))) stop("The length of x list, y list and rf model list should be identical."))
  perf_summ<-data.frame(matrix(NA, ncol=17, nrow=L*L)) 
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Accuracy", "AUC", "Kappa",
                         "Sensitivity", "Specificity", "Pos_Pred_Value","Neg_Pred_Value", "Precision", "Recall", 
                         "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")
  predicted<-matrix(list(), ncol=2, nrow=L*L)
  colnames(predicted)<-c("predicted", "probabilities")
  for(i in 1:L){
    y<-y_list[[i]]
    x<-x_list[[i]]
    try(if(nlevels(y)==1) stop("Less than one level in the subgroup for classification"))
    #rf_model<-randomForest(x, y, ntree=5000, importance=T)
    #oob<-rf.out.of.bag(x, y, nfolds=nfolds, verbose=verbose, ntree=ntree)
    oob<-rf_model_list[[i]]
    #---
    #  RF Training accuracy
    #---
    cat("\nTraining dataset: ", names(x_list)[i] ,"\n\n")
    conf<-confusionMatrix(data=oob$predicted, oob$y, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the self-validation: ", acc ,"\n") 
    #---
    #  AUC computation using "pROC" package
    #---
    auc<-get.auroc(oob$probabilities, oob$y, positive_class)
    cat("AUC in the self-validation: ", auc ,"\n") 
    
    D<-1:L
    T<-D[D!=i]
    a=1+(i-1)*L
    perf_summ[a, 1:3]<-c(names(x_list)[i], names(x_list)[i], "self_validation")
    perf_summ[a, 4:17]<-c(acc, auc, kappa_oob, conf$byClass)
    predicted[a, 1][[1]]<-data.frame(test_y=y, pred_y=oob$predicted)
    predicted[a, 2][[1]]<-oob$probabilities
    #rownames(predicted)[a]<-paste(names(x_list)[i], names(x_list)[i], sep="__VS__")
    loop_num<-1
    for(j in T){
      if(nrow(x_list[[j]])>0){
        newx<-x_list[[j]]
        newy<-y_list[[j]]
        #pred_prob<-predict(oob$rf.model, x_list[[j]], type="prob") # for regular rf.out.of.bag (randomForest) function
        pred_prob <- get.predict.probability.from.forest(oob$rf.model, newx) # ranger only
        pred_prob<-pred_prob[,order(colnames(pred_prob))] # to avoid unanticipated order of numeric levels of factor y
        pred_newy<-factor(predict(oob$rf.model, newx, type="response")$predictions) # ranger only
        if(identical(levels(newy), levels(oob$y))){
          colnames(pred_prob)<- levels(newy)
          levels(pred_newy)<- levels(newy)
          #---Accuracy
          cat("Test dataset: ", names(x_list)[j] ,"\n")
          test_conf<-confusionMatrix(data=pred_newy, newy, positive=positive_class)
          test_acc<-test_conf$overall[1]
          test_kappa<-test_conf$overall[2]
          cat("Accuracy in the cross-applications: ", test_acc ,"\n") 
          #---AUC
          test_auc<-get.auroc(pred_prob, newy, positive_class)
          cat("AUC in the cross-applications: ", test_auc ,"\n") 
          perf_summ[a+loop_num, 4:17]<-c(test_acc, test_auc, test_kappa, test_conf$byClass)
        }else{
          colnames(pred_prob)<- levels(oob$y)
          levels(pred_newy)<- levels(oob$y)
          perf_summ[a+loop_num, 4:17]<-rep(NA, 14)
        }
        perf_summ[a+loop_num, 1:3]<-c(names(x_list)[i], names(x_list)[j], "cross_application")
        predicted[a+loop_num, 1][[1]]<-data.frame(test_y=newy, pred_y=pred_newy)
        predicted[a+loop_num, 2][[1]]<-pred_prob
        #rownames(predicted)[a+loop_num]<-paste(names(x_list)[i], names(x_list)[j], sep="__VS__")
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  res
}

#' Based on pre-computed rf models regressing 'c_category' in each the sub-datasets splited by the 's_category',
#' perform cross-datasets application of the rf models. The inputs are precalculated 
#' rf regression models, and the outputs include accuracy, auc and Kappa statistics.
#' 
#' @param rf_model_list A list of rf.models generated from the function rf.out.of.bag.
#' @param df Training data: a data.frame.
#' @param metadata Sample metadata with at least two columns.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'numeric' sample label for rf regression in each of splited datasets.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_c=factor(c(rep("C", 15), rep("D", 15), rep("E", 15), rep("F", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_c', c_category='age')
#' rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
#' 

rf_reg.cross_appl<-function(rf_list, x_list, y_list){
  L<-length(rf_list$rf_model_list)
  try(if(!identical(L, length(x_list), length(y_list))) stop("The length of x list, y list and rf model list should be identical."))
  try(if(all(unlist(lapply(y_list, mode))!="numeric")) stop("All elements in the y list should be numeric for regression."))
  perf_summ<-data.frame(matrix(NA, ncol=14, nrow=L*L)) 
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Sample_size", "Min_acutal_value", "Max_acutal_value", "Min_predicted_value", "Max_predicted_value",
                         "MSE", "RMSE", "MAE", "MAE_perc", "R_squared", "Adj_R_squared")
  predicted<-list()
  for(i in 1:L){
    y<-y_list[[i]]
    x<-x_list[[i]]
    try(if(nlevels(y)==1) stop("Less than one level in the subgroup for classification"))
    oob<-rf_list$rf_model_list[[i]]
    #---  RF Training performance: MSE, MAE and R_squared
    cat("\nTraining dataset: ", names(x_list)[i] ,"\n\n")
    train_sample_size<-length(y)
    train_y_min<-range(y)[1] #paste0(range(y), collapse = "-")
    train_y_max<-range(y)[2]
    train_pred_y_min<-range(oob$predictions)[1]
    train_pred_y_max<-range(oob$predictions)[2]#paste0(range(oob$predictions), collapse = "-")
    train_MSE<-rf_list$rf_MSE[[i]]
    train_RMSE<-rf_list$rf_RMSE[[i]]
    train_MAE<-rf_list$rf_MAE[[i]]
    train_MAE_perc<-rf_list$rf_MAE_perc[[i]]
    train_R_squared<-rf_list$rf_R_squared[[i]]
    train_Adj_R_squared<-rf_list$rf_Adj_R_squared[[i]]
    cat("MSE in the self-validation: ", train_MSE ,"\n")
    cat("RMSE in the self-validation: ", train_RMSE ,"\n") 
    cat("MAE in the self-validation: ", train_MAE ,"\n") 
    cat("MAE percentage in the self-validation: ", train_MAE_perc ,"\n") 
    cat("R squared in the self-validation: ", train_R_squared ,"\n") 
    cat("Adjusted R squared in the self-validation: ", train_Adj_R_squared ,"\n") 
    D<-1:L
    T<-D[D!=i]
    a=1+(i-1)*L
    perf_summ[a, 1:3]<-c(names(x_list)[i], names(x_list)[i], "self_validation")
    perf_summ[a, 4:14]<-c(train_sample_size, train_y_min, train_y_max, train_pred_y_min, train_pred_y_max,
                          train_MSE, train_RMSE, train_MAE, train_MAE_perc, train_R_squared, train_Adj_R_squared)
    predicted[[a]]<-data.frame(test_y=y, pred_y=oob$predictions)
    names(predicted)[a]<-paste(names(x_list)[i], names(x_list)[i], sep="__VS__")
    loop_num<-1
    for(j in T){
      if(nrow(x_list[[j]])>0){
        newx<-x_list[[j]]
        newy<-y_list[[j]]
        pred_newy<-predict(rf_list$rf_model_list[[i]], newx, type="response")$predictions # ranger only
        #---  RF test performance: MSE, MAE and R_squared
        cat("Test dataset: ", names(x_list)[j] ,"\n")
        test_sample_size<-length(newy)
        test_y_min<-range(newy)[1] #paste0(range(newy), collapse = "-")
        test_y_max<-range(newy)[2]
        test_pred_y_min<-range(pred_newy)[1]
        test_pred_y_max<-range(pred_newy)[2]#paste0(range(pred_newy), collapse = "-")
        test_MSE<-MSE(newy, pred_newy) 
        test_RMSE<-RMSE(newy, pred_newy) 
        test_MAE<-MAE(newy, pred_newy) 
        test_MAE_perc<-MAE_perc(newy, pred_newy)
        test_R_squared<-R_squared(newy, pred_newy)
        test_Adj_R_squared<-Adj_R_squared(newy, pred_newy, k=ncol(newx))
        cat("MSE in the cross-applications: ", test_MSE ,"\n") 
        cat("RMSE in the cross-applications: ", test_RMSE ,"\n") 
        cat("MAE in the cross-applications: ", test_MAE ,"\n") 
        cat("MAE percentage in the cross-applications: ", test_MAE_perc ,"\n") 
        cat("R squared in the cross-application: ", test_R_squared ,"\n") 
        cat("Adjusted R squared in the cross-application: ", test_Adj_R_squared ,"\n") 
        perf_summ[a+loop_num, 1:3]<-c(names(x_list)[i], names(x_list)[j], "cross_application")
        perf_summ[a+loop_num, 4:14]<-c(test_sample_size, test_y_min, test_y_max, test_pred_y_min, test_pred_y_max,
                                       test_MSE, test_RMSE, test_MAE, test_MAE_perc, test_R_squared, test_Adj_R_squared)
        predicted[[a+loop_num]]<-data.frame(test_y=newy, pred_y=pred_newy)
        names(predicted)[a+loop_num]<-paste(names(x_list)[i], names(x_list)[j], sep="__VS__")
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  res
}

#' To generate sub-datasets and sub-metadata by pairing 
#' one level and each of other levels of the specified category in the metadata. 
#' 
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @return ...
#' 
generate.comps_datalist<-function(df, f, comp_group){
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  f_list<-list()
  df_list<-list()
  for(i in 1:L){
    f_list[[i]]<-factor(f[which(f==comp_group | f==all_other_groups[i])])
    df_list[[i]]<-df[which(f==comp_group | f==all_other_groups[i]), ]
  }
  names(f_list)<-names(df_list)<-paste(comp_group, all_other_groups, sep="_VS_")
  res<-list()
  res$f_list<-f_list
  res$df_list<-df_list
  res
}




#' Runs standard random forests with oob estimation for classification of 
#' one level VS all other levels of one category in the datasets. 
#' The output includes a list of rf models for the sub datasets
#' and all important statistics for each of features.
#' 
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. 
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))) 
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f, comp_group, verbose=FALSE, ntree=500, p.adj.method = "bonferroni", q_cutoff=0.05)                
#' comps_res
#' 
rf_clf.comps<-function(df, f, comp_group, verbose=FALSE, clr_transform=TRUE, rf_imp_values=FALSE, 
                       ntree=500, p.adj.method = "bonferroni", q_cutoff=0.05){
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  nCores <- detectCores()
  registerDoMC(nCores)
  oper<-foreach(i=1:L, .combine='comb', .multicombine=TRUE, .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {
    sub_f<-factor(f[which(f==comp_group | f==all_other_groups[i])])
    sub_df<-df[which(f==comp_group | f==all_other_groups[i]), ]
    print(levels(factor(sub_f)))
    # 1. sample size of all datasets
    sample_size<-length(factor(sub_f))
    dataset<-paste(comp_group, all_other_groups[i], sep="_VS_")
    # 2. AUC of random forest model
    oob <- rf.out.of.bag(sub_df, factor(sub_f), verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_values)
    rf_AUC <- get.auroc(oob$probabilities, factor(sub_f), comp_group)
    # 3. # of significantly differential abundant features between health and disease
    out<-BetweenGroup.test(sub_df, factor(sub_f), clr_transform=clr_transform, q_cutoff=q_cutoff, 
                           positive_class=comp_group, p.adj.method = p.adj.method)
    wilcox<-data.frame(feature=rownames(out),
                       dataset=rep(dataset, ncol(sub_df)), rf_imps=oob$importances, 
                       out)
    list(x=sub_df, y=sub_f, sample_size=sample_size, datasets=dataset, oob=oob, rf_AUC=rf_AUC, wilcox=wilcox)
  }
  names(oper[[1]])<-names(oper[[2]])<-names(oper[[5]])<-unlist(oper[[4]])
  result<-list()
  result$x_list<-oper[[1]]
  result$y_list<-oper[[2]]
  result$sample_size<-unlist(oper[[3]])
  result$datasets<-unlist(oper[[4]])
  result$rf_model_list<-oper[[5]]
  result$rf_AUC<-unlist(oper[[6]])
  result$feature_imps_list<-oper[[7]]
  
  return(result)
}

#' Runs standard random forests with oob estimation for classification of 
#' one level VS all other levels of one category in the datasets. 
#' The output includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' 
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. 
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param outdir The outputh directory, default is "./".
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))) 
#' comp_group="A"
#' rf_clf.comp.summ(df, f, comp_group, verbose=FALSE, ntree=500, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir="./")                
#' 
rf_clf.comp.summ<-function(df, f, comp_group, clr_transform=TRUE, nfolds=3, verbose=FALSE, ntree=5000, 
                           p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir="./"){
  res_list<-rf_clf.comps(df, f, comp_group, verbose=verbose, ntree=ntree, clr_transform=clr_transform)
  stopifnot(all(names(res_list) %in% c("datasets","rf_model_list","sample_size","rf_AUC","feature_imps_list")==TRUE))
  plot_res_list<-plot.res_list(res_list, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=outdir)
  result<-list()
  result$rf_models<-res_list$rf_model_list
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  return(result)
}

#' Non-specific features across datasets (at least present in two of datasets): Last update: 20190130
#' 
#' @param feature_res the inheritant output from the function of plot.res_list, rf_clf.comp.summ or rf_clf.by_dataset.summ 
#' @param level1 A string indicates the class of the factor, default is 'disease'
#' @param level2 A string indicates the class of the factor, default is 'health'
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param outdir The outputh directory, default is "./".
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' df <- data.frame(rbind(t(rmultinom(15, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299)))) # No feature significantly changed in all sub-datasets! 
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))) 
#' comp_group="A"
#' res <-rf_clf.comp.summ(df, f, comp_group, verbose=FALSE, ntree=500, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir="./")                
#' feature_res<-res$feature_res
#' id_non_spcf_markers(feature_res, positive_class="disease", other_class="health", p.adj.method= "BH", outdir="./")
#' 
id_non_spcf_markers <- function(feature_res, positive_class="disease", other_class="health", p.adj.method = "BH", outdir="./"){
  feature_res<-feature_res[order(feature_res$dataset, feature_res$feature), ]
  if(all(feature_res$Enr=='Neutral')) stop("No feature significantly changed in any sub-datasets!")
  Enr_dataset_df<-do.call(rbind, tapply(feature_res$Enr, feature_res$feature, summary))
  #Enr_dataset_df<-dcast(feature_res_sig, feature~Enr, fun=length, value.var="Enr")
  shared_global<-as.factor(apply(Enr_dataset_df, 1, function(x) ifelse(any(x[2]==sum(x) | x[3]==sum(x)), "all_shared", "unshared")))
  non_spcf_global<-as.factor(apply(Enr_dataset_df, 1, function(x) ifelse(any(x[2]==sum(x) | x[3]==sum(x)), "all_shared", 
                                                                         ifelse(any(x[2]>=2 | x[3]>=2), "non_spcf", 
                                                                                ifelse(any(x[3]==0 && x[2]==0), "non_marker", "spcf"))) 
  ))
  #---- statistics summary of non_spcf markers aross all datasets
  non_spcf_disease_global<-as.factor(apply(Enr_dataset_df, 1, function(x) 
    ifelse(any(x[3]>=2 && x[2]==0), paste("non_spcf_", positive_class, sep=""), 
           ifelse(any(x[3]==0 && x[2]>=2), paste("non_spcf_", other_class, sep=""), 
                  ifelse(any(x[3]>=1 && x[2]>=1), "non_spcf_mixed", 
                         ifelse(any(x[3]==0 && x[2]==0), "non_marker", "spcf"))))
  ))
  #---- calculate the overlap fraction of non-specific markers in each of datasets 
  tmp<-data.frame(feature_res, shared_global=as.character(shared_global), non_spcf_global=as.character(non_spcf_global), non_spcf_disease_global=as.character(non_spcf_disease_global), 
                  non_spcf_local=as.character(non_spcf_global), non_spcf_disease_local=as.character(non_spcf_disease_global), 
                  stringsAsFactors =F)
  tmp[which(tmp[, "Enr"]=="Neutral"),  c("non_spcf_local", "non_spcf_disease_local")]<-"non_marker"
  count_spcf<-do.call(rbind, tapply(tmp$non_spcf_local, tmp$dataset, function(x) 
    c(sum(x=="all_shared"), 
      sum(x=="non_spcf"),
      sum(x=="spcf"),
      sum(x!="all_shared" & x!="non_spcf" & x!="spcf")
    )
  )
  )
  zero_count <- tapply(tmp$mean_all, tmp$dataset, function(x) sum(x==0))
  count_spcf[, 4]<-count_spcf[, 4]-zero_count
  fac_spcf<-sweep(count_spcf, 1, rowSums(count_spcf), "/")
  colnames(fac_spcf)<-colnames(count_spcf)<-c("All_shared_markers","Non_specific_markers", "Specific_markers", "Others")
  # remove the columns with total zero values
  fac_spcf<-fac_spcf[, which(!apply(fac_spcf, 2, function(x) all(x==0)))]
  count_spcf<-count_spcf[, which(!apply(count_spcf, 2, function(x) all(x==0)))]
  fac_spcf<-data.frame(dataset=rownames(fac_spcf), fac_spcf)
  # output the fac_spcf table
  sink(paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".txt",sep=""));
  cat("\t"); write.table(fac_spcf,sep='\t',quote=F)
  sink(NULL)
  # output the count_spcf table
  sink(paste(outdir,"Markers_number_overlap_with_non_specific_",p.adj.method,".txt",sep=""));
  cat("\t"); write.table(count_spcf,sep='\t',quote=F)
  sink(NULL)
  # ggplot the fraction overlap of non-specific markers in all datasets
  fac_spcf_m<-melt(fac_spcf)
  fac_spcf_m$variable<-factor(fac_spcf_m$variable, levels = rev(levels(fac_spcf_m$variable)))
  mycolors<-c("grey80", gg_color_hue(nlevels(fac_spcf_m$variable)-1))
  p_summ<-ggplot(fac_spcf_m, aes(x=dataset, y=value, fill=variable)) + 
    geom_bar(stat="identity")+ 
    scale_fill_manual(values=mycolors)+
    ylab("Fraction overlap with non-specific markers")+
    coord_flip()+
    #ylim(0, 1)+ 
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".pdf",sep=""),plot=p_summ, width=6, height=4)
  #---- summary of the abundance and occurence rate of 
  #---- non-specific disease, non-specific health, non-specific mixed and non markers 
  #---- across all patients in all datasets
  feature_res_spcf<-tmp
  p_spcf_abd<-ggplot(feature_res_spcf, aes(x=non_spcf_disease_global, y=log(mean_all)))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("log10(mean abundance)")+
    coord_flip()+
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_specific_VS_mean_",p.adj.method,".pdf",sep=""),plot=p_spcf_abd, width=5, height=3)
  p_spcf_OccRate<-ggplot(feature_res_spcf, aes(x=non_spcf_disease_global, y=OccRate_all))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("Ubiquity")+
    coord_flip()+
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_specific_VS_ubiquity_",p.adj.method, ".pdf",sep=""),plot=p_spcf_OccRate, width=5, height=3)
  
  result<-list()
  result$fac_spcf<-fac_spcf
  result$count_spcf<-count_spcf
  result$feature_res_spcf<-feature_res_spcf
  result$plot_spcf_summ<-p_summ
  result$plot_spcf_abd<-p_spcf_abd
  result$plot_spcf_OccRate<-p_spcf_OccRate
  return(result)
  
}

