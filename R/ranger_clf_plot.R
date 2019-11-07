#--------------------------------------------------
p <- c("reshape2","randomForest", "optparse", "ade4", "doMC",
       "ggplot2", "RColorBrewer", "vegan", "xgboost", "caret")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))


#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y<-factor(c(rep("A", 20), rep("B", 10) ,rep("A", 10), rep("B", 20)))
#' system.time(rf.out.of.bag(x, y, imp_pvalues=FALSE))
#' system.time(rf.out.of.bag(x, y, imp_pvalues=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 1000))),
#'             t(rmultinom(8, 7500, rep(c(.201,.4,.12,.18,.099), 1000))),
#'             t(rmultinom(15, 7500, rep(c(.011,.3,.22,.18,.289), 1000))),
#'             t(rmultinom(15, 7500, rep(c(.091,.2,.32,.18,.209), 1000))),
#'             t(rmultinom(15, 7500, rep(c(.001,.1,.42,.18,.299), 1000)))))
#' y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' rf.out.of.bag(x_, y, imp_pvalues=FALSE)
#' #   user  system elapsed 
#' #100.824   0.263  42.947 
#' system.time(rf.out.of.bag(x, y, imp_pvalues=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7000, 75000, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(800, 7500, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(1500, 7500, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(1500, 7500, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(1500, 7500, c(.001,.1,.42,.18,.299)))))
#' y_ <- factor(c(rep("A", 1500), rep("B", 1500), rep("C", 1500), rep("D", 1500)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=TRUE))
#' y0<-factor(c(rep("old", 30), rep("young", 30)))
#' rf.out.of.bag(x, y0, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y0, imp_pvalues=TRUE)
#' y<- 1:60
#' rf.out.of.bag(x, y)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' 

plot_ROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir="./"){
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  pdf(paste(outpath, prefix, ".rf_clf_pROC.ci.pdf",sep=""), width=4, height=4)
  rocobj <- plot.roc(y, rf_clf_model$probabilities[, positive_class], main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
  ciobj <- ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
  plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
  text(70,25, paste0(levels(y), collapse = " VS "), pos=4)
  text(70,15, paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
  ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
  ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
  text(70,5, paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
  dev.off()
  result<-list()
  result$rocobj<-rocobj
  result$ciobj<-ciobj
  invisible(result)
}

plot_probabilities<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir="./"){
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  l<-levels(y); l_sorted<-sort(levels(y))
  Mycolor <- rep(c("#D55E00", "#0072B2"), length.out=length(l))
  if(identical(order(l), order(l_sorted))){
    Mycolor=Mycolor; l_ordered=l
  }else{Mycolor=rev(Mycolor); l_ordered=l_sorted}
  y_prob<-data.frame(y, rf_clf_model$probabilities)
  p<-ggplot(y_prob, aes(x=y, y=get(positive_class))) + 
    geom_violin()+
    geom_jitter(position=position_jitter(width=0.2),alpha=0.1) +
    geom_boxplot(outlier.shape = NA, width=0.4, alpha=0.01)+
    geom_hline(yintercept=0.5, linetype="dashed")+
    ylim(0, 1)+
    theme_bw()+
    #ggtitle(paste("Wilcoxon Rank Sum Test:\n P=",p_mf,sep=""))+
    xlab("") + 
    ylab(paste("Probability of ", positive_class))+
    theme(legend.position="none")+
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())+
    scale_color_manual(values = Mycolor, labels=l_ordered)
  ggsave(filename=paste(outpath,prefix,".probability_",positive_class,".boxplot.pdf",sep=""), plot=p, width=3, height=4)
}

plot_ClfPerf_VS_NumOfFeatures <- function(x, y, rf_clf_model, clf_perf="AUROC", positive_class=NA, outdir="./"){
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  if(rf_clf_model$error.type=="oob"){
    rf_imp_rank<-rank(-(rf_clf_model$importances)) #rf_all$importances
  }else{
    rf_imp_rank<-rank(-(rowMeans(rf_clf_model$importances))) #rf_all$importances
  }
  max_n<-max(rf_imp_rank, na.rm = TRUE)
  n_total_features<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[1:which.min(abs(n_features-n_total_features))]
  if(nlevels(y)==2){
    top_n_perf<-matrix(NA, ncol=5, nrow=length(n_features)+1)
    colnames(top_n_perf)<-c("n_features", "AUROC", "Accuracy", "Kappa", "F1")
    rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
    for(i in 1:length(n_features)){
      idx<-which(rf_imp_rank<=n_features[i])
      top_n_features<-names(rf_imp_rank[idx])
      x_n<-x[, top_n_features]
      y_n<-y
      top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500) # it depends on the rf.out.of.bag defined before
      top_n_AUC<-get.auroc(top_n_rf$probabilities, y_n, positive_class=positive_class)
      top_n_conf<-caret::confusionMatrix(data=top_n_rf$predicted, top_n_rf$y, positive=positive_class)
      top_n_perf[i, 1]<-n_features[i]
      top_n_perf[i, 2]<-top_n_AUC
      top_n_perf[i, 3]<-top_n_conf$overall[1] # Accuracy
      top_n_perf[i, 4]<-top_n_conf$overall[2]# kappa conf$byClass["F1"]
      top_n_perf[i, 5]<-top_n_conf$byClass["F1"]
    }
    all_AUC<-get.auroc(rf_clf_model$probabilities, y_n, positive_class=positive_class)
    all_conf<-caret::confusionMatrix(data=rf_clf_model$predicted, top_n_rf$y, positive=positive_class)
    top_n_perf[length(n_features)+1, ]<-c(max_n, all_AUC, all_conf$overall[1], 
                                          all_conf$overall[2], all_conf$byClass["F1"])
  }else{
    top_n_perf<-matrix(NA, ncol=3, nrow=length(n_features)+1)
    colnames(top_n_perf)<-c("n_features", "Accuracy", "Kappa")
    rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
    for(i in 1:length(n_features)){
      idx<-which(rf_imp_rank<=n_features[i])
      top_n_features<-names(rf_imp_rank[idx])
      x_n<-x[, idx]
      y_n<-y
      top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500) # it depends on the rf.out.of.bag defined before
      top_n_conf<-caret::confusionMatrix(data=top_n_rf$predicted, top_n_rf$y, positive=positive_class)
      top_n_perf[i, 1]<-n_features[i]
      top_n_perf[i, 2]<-top_n_conf$overall[1] # Accuracy
      top_n_perf[i, 3]<-top_n_conf$overall[2]# kappa 
    }
    all_conf<-caret::confusionMatrix(data=rf_clf_model$predicted, top_n_rf$y, positive=positive_class)
    top_n_perf[length(n_features)+1, ]<-c(max_n, all_conf$overall[1], all_conf$overall[2])
  }
  
  top_n_perf<-data.frame(top_n_perf)
  breaks<-top_n_perf$n_features
  p<-ggplot(top_n_perf, aes(x=n_features, y=get(clf_perf))) + 
    xlab("# of features used")+
    ylab(clf_perf)+
    scale_x_continuous(trans = "log",breaks=breaks)+
    geom_point() + geom_line()+ ylim(0.5, 1)+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"), 
          panel.border = element_blank())
  ggsave(filename=paste(outpath,"train.",clf_perf,"_VS_top_ranking_features.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
}


