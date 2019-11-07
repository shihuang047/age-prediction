#' Loading libraries
p <- c("ade4","caret","reshape2","MASS","ranger","plyr", "fmsb", "foreach", "doMC", "ALEPlot", "vegan", "compositions",
       "microbenchmark", "lattice","gplots","ggplot2","squash","RColorBrewer","pheatmap","pROC","ROCR","caTools")
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
    suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))


#' runs cross-validation 
#' if predict.fun is NULL, uses S3 predict method
#' if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
#' ...: additional parameters for train.fun
#'
#' value:
#' y: true values
#' predicted: cv predicted values
#' probabilities: cv predicted class probabilities (or NULL if unavailable)
#' confusion.matrix: confusion matrix (true x predicted)
#' nfolds: nfolds
#' params: list of additional parameters
#' importances: importances of features as predictors

#' Get probability of mislabeling by several measures
#' returns matrix of p(alleged), max(p(others)), p(alleged) - max(p(others))

#rf.opts<-list(outdir="./", ntree=5000, verbose=FALSE, nfolds=3)

"get.mislabel.scores" <- function(y, y.prob){
    result <- matrix(0,nrow=length(y),ncol=3)
    # get matrices containing only p(other classes), and containing only p(class)
    mm <- model.matrix(~0 + y)
    y.prob.other.max <- apply(y.prob * (1-mm),1,max)
    y.prob.alleged <- apply(y.prob * mm, 1, max)
    result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
    rownames(result) <- rownames(y.prob)
    colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
    return(result)
}

#' Get balanced folds where each fold has close to overall class ratio
#' @param y Data label for each row.
#' @param nfolds The number of folds in the cross validation. 
#' @return folds
#' @examples
#' 
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' print(balanced.folds(y, nfolds=3))
#' y<- 1:60
#' print(balanced.folds(y, nfolds=3))
#' 
"balanced.folds" <- function(y, nfolds=3){
    folds = rep(0, length(y))
    if(is.factor(y)){
      classes<-levels(y)
    }else{
      y<-factor(findInterval(y, quantile(y, seq(0, 1, by=1/nfolds), type=5), rightmost.closed=TRUE))
      classes<-levels(y)
    } 
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }
    else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}

#' Runs standard random forests with n-folds cross-validation error estimation
#' This is merely a wrapper that extracts relevant info
#'
#' @param x Training data: data.matrix or data.frame.
#' @param y Data label for each row.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. 
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' 
#' x0 <- data.frame(t(rmultinom(16,160,c(.001,.5,.3,.3,.299))) + 0.65) 
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y0<-factor(c(rep("A", 10), rep("B", 30), rep("C", 5), rep("D", 15)))
#' system.time(rf.cross.validation(x, y, imp_pvalues=FALSE))
#' system.time(rf.cross.validation(x, y, imp_pvalues=TRUE))
#' rf.cross.validation(x, y0, imp_pvalues=FALSE)
#' y<- 1:60
#' rf.cross.validation(x, y, nfolds=5, imp_pvalues=FALSE)
#' rf.cross.validation(x, y, nfolds=5, imp_pvalues=TRUE)

"rf.cross.validation" <- function(x, y, nfolds=3, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y, nfolds=nfolds)
    result <- list()
    if(is.factor(y)){
      result$y <- as.factor(y)
      result$errs <- numeric(length(unique(folds)))
      result$probabilities <- matrix(0, nrow=length(result$y), ncol=length(levels(result$y)))
      rownames(result$probabilities) <- rownames(x)
      colnames(result$probabilities) <- levels(result$y)
    }else{
      result$y <- y
      result$MSE <- numeric(length(unique(folds)))
      result$MAE <- numeric(length(unique(folds)))
      result$R_squared <- numeric(length(unique(folds)))
    }
    result$models<-list()
    result$predicted <- result$y
    result$importances <- matrix(0, nrow=ncol(x), ncol=nfolds)
    if(imp_pvalues==TRUE){ result$importance_pvalues <- matrix(0, nrow=ncol(x), ncol=nfolds) }
    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        cat("The # of folds: ", fold, "\n") # if(verbose)
        foldix <- which(folds==fold)
        if(is.factor(y)) y_tr<-factor(result$y[-foldix]) else y_tr<-result$y[-foldix]
        data<-data.frame(y=y_tr, x[-foldix,])
        require(Matrix)
        if(sparse){
          sparse_data <- Matrix(data.matrix(data), sparse = TRUE)
          result$models[[fold]]<- model <- ranger::ranger(dependent.variable.name="y", data=sparse_data, classification=ifelse(is.factor(y_tr), TRUE, FALSE), 
                                                          keep.inbag=TRUE, importance='permutation', verbose=verbose, num.trees=ntree)
        }else{
          result$models[[fold]]<- model <- ranger::ranger(y~., data=data, keep.inbag=TRUE, importance='permutation', 
                                                         classification=ifelse(is.factor(y_tr), TRUE, FALSE), num.trees=ntree, verbose=verbose)
        }
        
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx, nrow=1)
        predicted_foldix<-predict(model, newx)$predictions
        if(is.factor(y)){
          if(sparse){
            y_numeric<-as.factor(sparse_data[,'y'])
            predicted_foldix <- factor(predicted_foldix, levels=levels(y_numeric))
            levels(predicted_foldix)=levels(y)
          }else{
            predicted_foldix <- factor(predicted_foldix, levels=levels(y))
          }
          result$predicted[foldix] <- predicted_foldix
          probs <- get.predict.probability.from.forest(model, newx); colnames(probs)<-levels(result$y)
          result$probabilities[foldix, colnames(probs)] <- probs
          result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
          result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
          cat("Error rate: ", result$errs[fold], "\n")
        }else{
          result$predicted[foldix] <- predicted_foldix
          reg_perf<-get.reg.predict.performance(model, newx, newy=result$y[foldix])
          result$MSE[fold] <- reg_perf$MSE 
          result$RMSE[fold] <- reg_perf$RMSE
          result$MAE[fold] <- reg_perf$MAE 
          result$MAE_perc[fold] <- reg_perf$MAE_perc
          result$R_squared[fold] <- reg_perf$R_squared 
          result$Adj_R_squared[fold] <- reg_perf$Adj_R_squared 
          cat("Mean squared residuals: ", reg_perf$MSE, "\n")
          cat("Mean absolute error: ", reg_perf$MAE, "\n")
          cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
        }
        if(imp_pvalues==FALSE){
          result$importances[,fold] <- model$variable.importance
          }else{
          imp<-importance_pvalues(model, method = "altmann", formula = y ~., data=data)
          result$importances[,fold] <- imp[, 1]
          result$importance_pvalues[,fold] <- imp[, 2]
        }
        
    }
    result$nfolds <- nfolds
    result$params <- list(...)
    result$error.type <- "cv"
    return(result)
}

#' Runs standard random forests with out-of-bag error estimation
#' This is merely a wrapper that extracts relevant info
#' Return values are the same as rf.cross.validation
#'
#' @param x Training data: data.matrix or data.frame.
#' @param y Data label for each row.
#' @param ntree The number of trees. 
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return ...
#'
#' @seealso ranger
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

"rf.out.of.bag" <-function(x, y, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
    set.seed(123)
    data<-data.frame(y, x)
    require(Matrix)
    if(sparse){
      sparse_data <- Matrix(data.matrix(data), sparse = sparse)
      rf.model <- ranger::ranger(dependent.variable.name="y", data=sparse_data, keep.inbag=TRUE, importance='permutation', 
                                 classification=ifelse(is.factor(y), TRUE, FALSE), num.trees=ntree, verbose=verbose, probability = FALSE)
    }else{
      rf.model<-ranger::ranger(y~., data=data, keep.inbag=TRUE, importance='permutation', 
                               classification=ifelse(is.factor(y), TRUE, FALSE), num.trees=ntree, verbose=verbose, probability = FALSE)
    }
    result <- list()
    result$rf.model <- rf.model
    result$y <- y
    if(is.factor(y)){
      if(sparse){
        y_numeric<-as.factor(sparse_data[,'y'])
        result$predicted <- factor(rf.model$predictions,levels=levels(y_numeric)); levels(result$predicted)<-levels(y)
      }else{
        result$predicted <- factor(rf.model$predictions,levels=levels(y))
      }
      result$probabilities <- get.oob.probability.from.forest(rf.model, x); colnames(result$probabilities)<-levels(result$y)
      result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
      result$errs <- mean(result$predicted != result$y)
      #if(nlevels(y)==2) result$auroc <- get.rf.auroc(result$probabilities, y, positive_class = levels(y)[1])
    }else{
      result$predicted <- rf.model$predictions
      reg_perf<-get.reg.oob.performance(rf.model, y)
      result$MSE <- reg_perf$MSE
      result$RMSE <- reg_perf$RMSE 
      result$MAE <- reg_perf$MAE 
      result$MAE_perc <- reg_perf$MAE_perc 
      result$R_squared <- reg_perf$R_squared 
      result$Adj_R_squared <- reg_perf$Adj_R_squared
      cat("Mean squared residuals: ", reg_perf$MSE, "\n")
      cat("Mean absolute error: ", reg_perf$MAE, "\n")
      cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
    }
    result$params <- list(ntree=rf.opts$ntree)
    if(imp_pvalues==FALSE){
      result$importances <- rf.model$variable.importance
      result$importances[colSums(x)==0]<-NA
      }else{
      result$importances <- importance_pvalues(rf.model, method = "altmann", formula = y ~., data=data)
      }
    result$error.type <- "oob"
    return(result)
}

#' get probability of each class using only out-of-bag predictions from RF
#' @param model A object of random forest out-of-bag model.
#' @param x The training data.
#' @return probabilities
#' 
"get.oob.probability.from.forest" <- function(model, x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- get.oob.votes.from.forest(model,x)
    # convert to probs
    probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
    rownames(probs) <- rownames(x)
    colnames(probs) <- model$forest$class.values
    
    return(invisible(probs))
}

#' get votes for each class using only out-of-bag predictions from RF
#' @param model A object of random forest out-of-bag model.
#' @param x The training data.
#' @return votes
#' 
"get.oob.votes.from.forest" <- function(model, x){
    # get aggregated class votes for each sample using only OOB trees
    n_classes <- length(model$forest$class.values)
    votes <- matrix(0, nrow=nrow(x), ncol=n_classes)
    
    rf.pred <- predict(model, data.frame(x), type="response", predict.all=TRUE)
    model_inbag<-t(do.call(rbind, model$inbag)) # convert the inbag list (length=ntree) into a inbag matrix
    for(i in 1:nrow(x)){
        # find which trees are not inbag for this sample
        outofbag <- model_inbag[i,]==0
        # get oob predictions for this sample
        votes[i,] <- table(factor(rf.pred$predictions[i,][outofbag],levels=seq_along(model$forest$class.values)))
    }
    rownames(votes) <- rownames(x)
    colnames(votes) <- model$forest$levels
    
    return(invisible(votes))
}
#' to get votes from the forest in the external test. Last update: Apr. 7, 2019
#' @param model A object of random forest out-of-bag model.
#' @param newx The external test data.
#' @return probabilities
#' 
"get.predict.probability.from.forest" <- function(model, newx){
  # get aggregated class votes for each sample using only OOB trees
  n_classes <- length(model$forest$class.values) #model$forest$levels
  votes <- matrix(0, nrow=nrow(newx), ncol=n_classes)
  rf.pred <- predict(model, data.frame(newx), type="response", predict.all=TRUE)
  for(i in 1:nrow(newx)){
    votes[i,] <- table(factor(rf.pred$predictions[i,], levels=seq_along(model$forest$class.values)))
  }
  rownames(votes) <- rownames(newx)
  colnames(votes) <- seq_along(model$forest$class.values) # model$forest$levels; to correct unanticipated order of class values
  probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
  return(invisible(probs))
}

#' to get oob performance from forests for training data. Last update: Apr. 7, 2019
#' @param model A object of random forest out-of-bag regression model.
#' @param y The dependent variable in training data.
#' @return perf
#' 

"get.reg.oob.performance" <- function(model, y){
  pred_y<-as.numeric(model$predictions)
  y<-as.numeric(y)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(mean((y-pred_y)^2)) 
  MAE <- mean(sqrt((y-pred_y)^2))
  MAE_perc <- mean(sqrt((y-pred_y)^2)/y)
  R_squared <- 1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))
  Adj_R_squared <- Adj_R_squared(y, pred_y, k = model$num.independent.variables)
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$MAE<-MAE
  perf$MAE_perc<-MAE_perc
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}

#' Regression performance metrics
#' @example 
#' y = c(1, 2, 4, 8)
#' pred_y = c(1.1, 1.8, 4.5, 8)
#' pred_y = c(14, 8, 5, 2)
#' MSE(y, pred_y)
#' RMSE(y, pred_y)
#' MAE(y, pred_y)
#' MAE_perc(y, pred_y)
#' R_squared(y, pred_y)
#' Adj_R_squared(y, pred_y, 200)

MSE<-function(y, pred_y){ mean((y-pred_y)^2)}
RMSE<-function(y, pred_y){ sqrt(mean((y-pred_y)^2))}
MAE<-function(y, pred_y){ mean(sqrt((y-pred_y)^2))}
R_squared<-function(y, pred_y){ 1-(sum((y-pred_y)^2) / sum((y-mean(y))^2)) }
MAE_perc<-function(y, pred_y){ mean(sqrt((y-pred_y)^2)/y)}
Adj_R_squared<-function(y, pred_y, k){
  n=length(y); 
  1-(1-R_squared(y, pred_y)^2)*(n-1)/(n-k-1) # k is # of predictors
} 


#' to get prediction performance from forests for the external test data. Last update: Apr. 7, 2019
#' @param model A object of random forest out-of-bag regression model.
#' @param newx The external test data.
#' @param newy The dependent variable in the external test data.
#' @return perf
#' 
"get.reg.predict.performance" <- function(model, newx, newy){
  rf.pred <- predict(model, data.frame(newx))
  pred_y <- as.numeric(rf.pred$predictions)
  y <- as.numeric(newy)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(mean((y-pred_y)^2)) 
  MAE <- mean(sqrt((y-pred_y)^2))
  MAE_perc <- mean(sqrt((y-pred_y)^2)/y)
  R_squared <- 1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))
  Adj_R_squared <- Adj_R_squared(y, pred_y, k = ncol(newx))
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$MAE<-MAE
  perf$MAE_perc<-MAE_perc
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}

#' 'get.auroc' calculates the area under the ROC curve
#' @param obs A binary factor vector indicates observed values.
#' @param prob A two-column numeric matrix of the same # of rows as the length of observed values,
#'             containing the predicted value of each observation. 
#' @return auroc
#' @example 
#' 
#' obs<-factor(c(rep("A", 31), rep("B", 29)))
#' pred <- c(runif(30, 0.5, 0.9), runif(30, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' get.auroc(obs, prob, positive_class="A")
#' 
get.auroc <- function(prob, obs, positive_class) {
  pos_class<-function(f, positive_class){
    if(nlevels(f)>2) stop("Only two-level factor allowed.")
    idx<-which(levels(f)==positive_class)
    levels(f)[-idx]<-0
    levels(f)[idx]<-1
    factor(f)
  }
  obs<-pos_class(obs, positive_class=positive_class)
  require(ROCR)
  pred <- prediction(prob[, positive_class], obs)
  auroc  <- performance(pred, "auc")@y.values[[1]]
  auroc
  return(auroc)
}

#' 'get.auprc' calculates the area under Precision-recall curve
#' @param obs A binary factor vector indicates observed values.
#' @param prob A two-column numeric matrix of the same # of rows as the length of observed values,
#'             containing the predicted value of each observation. 
#' @example 
#' 
#' obs<-factor(c(rep("A", 10), rep("B", 50)))
#' pred <- c(runif(10, 0.4, 0.9), runif(50, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' get.auprc(obs, prob, "A")
#' get.auprc(obs, prob, "B")
#' 
get.auprc <- function(obs, prob, positive_class) {
  pos_class<-function(f, positive_class){
    if(nlevels(f)>2) stop("Only two-level factor allowed.")
    idx<-which(levels(f)==positive_class)
    levels(f)[-idx]<-0
    levels(f)[idx]<-1
    factor(f)
  }
  obs<-pos_class(obs, positive_class=positive_class)
  require(ROCR)
  require(caTools)
  xx.df <- prediction(prob[, positive_class], obs)
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  
  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))
  
  # Designate recall = 0 as precision = x...arbitrary
  xy <- rbind(c(0, 0), xy)
  #xy <- xy[!(rowSums(xy)==0), ]
  
  res   <- caTools::trapz(xy$recall, xy$precision)
  res
}

#' prints random forests results file
"save.rf_clf.results" <- function(result, feature.ids, rf.opts){
    save.rf_clf.results.summary(result, rf.opts, outdir=rf.opts$outdir)
    save.rf_clf.results.probabilities(result, outdir=rf.opts$outdir)
    save.rf_clf.results.mislabeling(result, outdir=rf.opts$outdir)
    save.rf_clf.results.importances(result, feature.ids, outdir=rf.opts$outdir)
    save.rf_clf.results.confusion.matrix(result, outdir=rf.opts$outdir)
}

#' Print "summary" file
"save.rf_clf.results.summary" <- function(result, rf.opts, filename='train.summary.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    err <- mean(result$errs)
    err.sd <- sd(result$errs)
    y<-result$y
    baseline.err <- 1-max(table(y))/length(y)
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat(sprintf('Model\tRandom Forest\n'))
    cat(sprintf('Error type\t%s\n',result$error.type))
    if(result$error.type == 'oob'){
        cat(sprintf('Estimated error\t%.5f\n',err))
    } else {
        cat(sprintf('Estimated error (mean +/- s.d.)\t%.5f +/- %.5f\n',err,err.sd))
    }
    cat(sprintf('Baseline error (for random guessing)\t%.5f\n',baseline.err))
    cat(sprintf('Ratio baseline error to observed error\t%.5f\n',baseline.err / err))
    cat(sprintf('Number of trees\t%d\n',result$params$ntree))
    sink(NULL)
}

#' Print "probabilities" file
"save.rf_clf.results.probabilities" <- function(result, filename='train.cv_probabilities.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(result$probabilities,sep='\t',quote=F)
    sink(NULL)
}

#' Print "mislabeling" file
"save.rf_clf.results.mislabeling" <- function(result, filename='train.mislabeling.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(get.mislabel.scores(result$y,result$probabilities),sep='\t',quote=F)
    sink(NULL)
}

"save.rf_clf.results.importances" <- function(result, feature.ids, filename='train.feature_importance_scores.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    if(is.null(dim(result$importances))){
        imp <- result$importances
        imp.sd <- rep(NA,length(imp))
    } else {
        imp <- rowMeans(result$importances)
        imp.sd <- apply(result$importances, 1, sd)
    }
    output.table <- cbind(imp, imp.sd)
    rownames(output.table) <- feature.ids
    output.table <- output.table[sort(imp,dec=T,index=T)$ix,]
    colnames(output.table) <- c('Mean_decrease_in_accuracy','Standard_deviation')

    sink(filepath)
    cat('Feature_id\t')
    write.table(output.table,sep='\t',quote=F)
    sink(NULL)
}

#' Print "confusion matrix" file
"save.rf_clf.results.confusion.matrix" <- function(result, filename='train.confusion_matrix.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    # add class error column to each row
    x <- result$confusion.matrix
    class.errors <- rowSums(x * (1-diag(nrow(x)))) / rowSums(x)
    output <- cbind(result$confusion.matrix, class.errors)
    colnames(output)[ncol(output)] <- "Class error"
    sink(filepath)
    cat('True\\Predicted\t')
    write.table(output,quote=F,sep='\t')
    sink(NULL)
}
