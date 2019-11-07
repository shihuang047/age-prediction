#--------------------------------------------------
p <- c("reshape2","randomForest", "optparse", "ade4", "biom", "doMC",
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
#' system.time(rf.out.of.bag(x, y, imp_pvalues=FALSE))
#' system.time(rf.out.of.bag(x, y, imp_pvalues=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 10000))),
#'             t(rmultinom(8, 7500, rep(c(.201,.4,.12,.18,.099), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.011,.3,.22,.18,.289), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.091,.2,.32,.18,.209), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.001,.1,.42,.18,.299), 10000)))))
#' y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' #   user  system elapsed 
#' #100.824   0.263  42.947 
#' system.time(rf.out.of.bag(x, y, imp_pvalues=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7000, 75000, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8000, 75000, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15000, 75000, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15000, 75000, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15000, 75000, c(.001,.1,.42,.18,.299)))))
#' y_ <- factor(c(rep("A", 15000), rep("B", 15000), rep("C", 15000), rep("D", 15000)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=TRUE))
#' y0<-factor(c(rep("old", 30), rep("young", 30)))
#' rf.out.of.bag(x, y0, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y0, imp_pvalues=TRUE)
#' y<- 1:60
#' rf.out.of.bag(x, y)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' 

plot_obs_VS_pred <- function(y, predicted_y, prefix="train", target_field, span=1, outdir="./"){
  df<-data.frame(y, predicted_y)
  p<-ggplot(df, aes(x=y, y=predicted_y))+
    ylab(paste("Predicted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_smooth(method="loess",span=span)+
    theme_bw()
  #coord_flip()+
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  sink(paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""));cat("\t");write.table(df, quote=FALSE,sep="\t");sink()
  invisible(p)
}

plot_residuals <- function(y, predicted_y, prefix="train", target_field, outdir="./"){
  df<-data.frame(y, predicted_y, rsdl=y-predicted_y)
  p<-ggplot(df, aes(x=y, y=rsdl))+
    ylab(paste("Residuals of prediceted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_hline(yintercept=0)+
    #geom_smooth(method="loess",span=span)+
    theme_bw()
  #coord_flip()+
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_residuals_of_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  sink(paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""));cat("\t");write.table(df, quote=FALSE,sep="\t");sink()
  invisible(p)
}


rand_MAE <- function(y, permutation=100){
  set.seed(123)
  rand_y_mat <-replicate(permutation, sample(y, replace = FALSE))
  apply(rand_y_mat, 2, function(x) MAE(y, x))
}

plot_MAE_VS_rand<-function(y, predicted_y, prefix="train", target_field, permutation, outdir="./"){
  MAE_value<-MAE(y, predicted_y)
  MAE_values<-data.frame(MAE_value, rand_MAE_values=rand_MAE(y))
  p<-ggplot(MAE_values, aes(x=rand_MAE_values)) + geom_histogram(alpha=0.5) + 
    #xlim(c(1, 20))+ 
    xlab("MAE")+
    geom_vline(data=MAE_values, aes(xintercept = MAE_value)) +
    annotate(geom="text", x=MAE_value, y=20, label=as.character(round(MAE_value, 2)), color="red", hjust = 0)
  theme_bw()
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".MAE_vs_rand.histogram.pdf",sep=""), plot=p, height=4, width=4)
  invisible(p)
}

plot_train_vs_test<-function(train_y, predicted_train_y, test_y, predicted_test_y, 
                             train_prefix="train", test_prefix="test", train_target_field, test_target_field, outdir="./"){
  train.pred<-data.frame(value=train_y,predicted_value=predicted_train_y)
  test.pred<-data.frame(value=test_y,predicted_value=predicted_test_y)
  data_name<-c(rep(train_prefix,nrow(train.pred)),rep(test_prefix,nrow(test.pred)))
  pred<-data.frame(data=data_name,rbind(train.pred,test.pred))
  pred$data<-factor(pred$data,levels=c(train_prefix,test_prefix),ordered=TRUE)
  l<-levels(pred$data); l_sorted<-sort(levels(pred$data))
  Mycolor <- c("#D55E00", "#0072B2") 
  if(identical(order(l), order(l_sorted))){Mycolor=Mycolor }else{Mycolor=rev(Mycolor)}
  p<-ggplot(pred,aes(x=value,y=predicted_value))+
    ylab(paste("Predicted ",train_target_field,sep=""))+
    xlab(paste("Observed ",train_target_field,sep=""))+
    geom_point(aes(color=data), alpha=0.1)+
    geom_smooth(aes(color=data), method="loess",span=1)+
    scale_color_manual(values = Mycolor)+
    theme_bw() +
    facet_wrap(~data)+
    theme(legend.position="none")
  #coord_flip()+
  ggsave(filename=paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_ggplot.pdf",sep=""),plot=p, height=4, width=8)
  sink(paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_results.xls",sep=""));
  cat("\t");write.table(pred, quote=FALSE,sep="\t");sink()
  invisible(p)
}

plot_RegPerf_VS_NumOfFeatures <- function(x, y, rf_model, reg_perf="MAE", outdir="./"){
  rf_imp_rank<-rank(-(rf_model$importances)) #rf_all$importances
  max_n<-max(rf_imp_rank, na.rm = TRUE)
  n_total_features<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[1:which.min(abs(n_features-n_total_features))]
  top_n_perf<-matrix(NA, ncol=5, nrow=length(n_features)+1)
  colnames(top_n_perf)<-c("n_features", "MSE", "RMSE", "MAE", "MAE_perc")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  for(i in 1:length(n_features)){
    idx<-which(rf_imp_rank<=n_features[i])
    top_n_features<-names(rf_imp_rank[idx])
    x_n<-x[, top_n_features]
    y_n<-y
    top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500) # it depends on the rf.out.of.bag defined before
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-top_n_rf$MSE
    top_n_perf[i, 3]<-top_n_rf$RMSE
    top_n_perf[i, 4]<-top_n_rf$MAE
    top_n_perf[i, 5]<-top_n_rf$MAE_perc
  }
  top_n_perf[11, ]<-c(max_n, rf_model$MSE, rf_model$RMSE, rf_model$MAE, rf_model$MAE_perc)
  top_n_perf<-data.frame(top_n_perf)
  breaks<-top_n_perf$n_features
  p<-ggplot(top_n_perf, aes(x=n_features, y=get(reg_perf))) + 
    xlab("# of features used")+
    ylab(paste(reg_perf, " (yrs)", sep=""))+
    scale_x_continuous(trans = "log",breaks=breaks)+
    geom_point() + geom_line()+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"), 
          panel.border = element_blank())
  ggsave(filename=paste(outpath,"rf__",reg_perf,"__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
}

calc_rel_predicted<-function(train_y, predicted_train_y, test_y=NULL, predicted_test_y=NULL, 
                             train_prefix="train", test_prefix="test", 
                             train_target_field="y", test_target_field=NULL, outdir="./"){
  spl_train <- smooth.spline(train_y, predicted_train_y)
  train_relTrain <- residuals(spl_train); names(train_relTrain)<-names(train_y)
  relTrain_data<-train_relTrain_data <- data.frame(y=train_y, predicted_y=predicted_train_y, rel_predicted_y=train_relTrain)
    sink(paste(outdir, train_prefix,".Relative_",train_target_field,".results.xls",sep=""));cat("\t");write.table(relTrain_data,quote=FALSE,sep="\t");sink()

  if(!is.null(test_y) & !is.null(predicted_test_y)){
    test_relTrain <- predicted_test_y - predict(spl_train, test_y)$y
    test_relTrain_data <- data.frame(y=test_y, predicted_y=predicted_test_y, rel_predicted_y=test_relTrain)
    relTrain_data<-rbind(train_relTrain_data, test_relTrain_data)
    DataSet <- factor(c(rep(train_prefix,length(train_relTrain)), rep(test_prefix,length(test_relTrain)) ))
    relTrain_data<-data.frame(relTrain_data, DataSet)
    sink(paste(outdir, train_prefix,"-",test_prefix,".Relative_",train_target_field,".results.xls",sep=""));
    cat("\t");
    write.table(relTrain_data,quote=FALSE,sep="\t");sink()
    
  }
    return(relTrain_data)
}

plot_rel_predicted <- function(relTrain_data, DataSet=NA, prefix="train", target_field, outdir="./"){
  if(!is.na(DataSet)){
    relTrain_data<-subset(relTrain_data, DataSet==DataSet)
    prefix<-DataSet
    p<-ggplot(relTrain_data, aes(x=y, y=rel_predicted_y))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      #geom_smooth(method="loess",span=span)+
      theme_bw()
    #coord_flip()+
  }else{
    p<-ggplot(relTrain_data, aes(x=y, y=rel_predicted_y, color=DataSet))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      #geom_smooth(method="loess",span=span)+
      theme_bw()
    prefix="train-test"
  }
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_relative_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  invisible(p)
}


boxplot_rel_predicted_train_vs_test<-function(relTrain_data, train_target_field="y", 
                                              train_prefix="train", test_prefix="test", outdir="./"){
    NoTestDataset<-all(grepl("DataSet", colnames(relTrain_data))==FALSE)
    if(NoTestDataset) stop("Test dataset should be included for residuals comparison between train and test datasets!")
    p_w<-formatC(wilcox.test(rel_predicted_y~DataSet, data = relTrain_data)$p.value,digits=4,format="g")
    # l<-levels(relTrain_data$DataSet); l_sorted<-sort(levels(relTrain_data$DataSet))
    # Mycolor <- rep(c("#D55E00", "#0072B2"), length.out=length(l))
    # if(identical(order(l), order(l_sorted))){
    #   Mycolor=Mycolor; l_ordered=l
    # }else{Mycolor=rev(Mycolor); l_ordered=l_sorted}
    p<-ggplot(relTrain_data, aes(x=DataSet, y=rel_predicted_y)) + 
      geom_violin(aes(color=DataSet))+
      geom_boxplot(outlier.shape = NA, width=0.4)+
      geom_jitter(position=position_jitter(width=0.2),alpha=0.1) + # aes(color=DataSet),
      geom_hline(yintercept=0)+
      theme_bw()+
      ggtitle(paste("Wilcoxon Rank Sum Test:\n P=", p_w, sep=""))+
      xlab("Data sets") + 
      ylab(paste("Relative microbiota", train_target_field))+
      theme(legend.position="none")
    # p<-p+scale_color_manual(values = Mycolor, labels=l_ordered)
    ggsave(filename= paste(outdir, train_prefix,"-",test_prefix, ".Relative_",train_target_field,".boxplot.pdf",sep=""), width=3, height=4)
    invisible(p)
}

