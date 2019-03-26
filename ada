library("rpart")
library("rpart.plot")
library("crayon")
####################################################################################################
####################################################################################################
data <- read.csv("audit_risk.csv", TRUE)
X_where = 1:5
y_where = 28
how_many = 20
cm_res <- list()
bm_res <- list()
xx <- 1:how_many
X <- data[, X_where]
y <- data[, y_where]
####################################################################################################
####################################################################################################

adaboost <- function(X, y, cm, n_rounds = 100, 
                     control = rpart.control(cp = -1, maxdepth = 1)){
  
  # count the number of rows
  n      <- nrow(X)  
  
  
  # assign the weight
  if(cm == TRUE){
    w <- rep(1/n, n)
    #w <- rep(1/n, n) * sqrt(exp(225/306)) ,, 126/306
    w <- ifelse(y == 1, (1/n)*(sqrt(exp(length((which(y == 1)))/nrow(X)))), (1/n)/(sqrt(exp((length(which(y == -1)))/nrow(X)))))
    #w <- ifelse(y == 1, (1/n)*(sqrt(exp(length((which(y == 1)))/nrow(X)))), (1/n))
    #w <- ifelse(y == 1, (1/n)*((exp(length((which(y == 1)))/nrow(X)))), (1/n)/((exp((length(which(y == -1)))/nrow(X)))))
    
    
    
  }
  
  else{
    
    w      <- rep(1/n, n)
    
  }
  
  # trees constructed in each round of boosting
  trees  <- list()
  
  # learning rate(weight) of each classifiers
  alphas <- list()
  
  rules <- list()
  
  if(cm == TRUE){
    for(i in seq(n_rounds)){
      
      tree <- rpart::rpart(y ~ ., 
                           data = data.frame(X), weights = w,
                           method = "class", control = control, 
                           x = FALSE, y = FALSE, model = TRUE)
      
      rule <- rpart.rules(tree)
      pred <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
      
      # weighted error rate
      e    <- sum(w * (pred != y))
      
      if(e >= 0.5){
        
        e <- 1 - e
      }
      
      else{
        
        e <- e
        
      }
      # learning rate(weight) of each classifiers
      alpha <- 1/2 * log((1-e)/e)
      #alpha <- log(-(1- ((1-e)/e)))
      #alpha <- (1/2)*log(((1-e)/e),50)
      
      # update weight of each data
      w     <- w * (exp((-alpha*pred*y)))
      #w     <- w  * ((1/(1+exp(alpha*pred*y)))*2)
      #w <- w * (50^((-alpha*pred*y)))
      # normalize weight of each data
      w     <- w / sum(w)
      
      # If classifier's error rate is nearly 0, boosting process ends
      if(abs(e) < 1e-5){
        
        # if first base classifier predicts data perfectly, boosting process ends
        if(i == 1){
          
          trees[[i]]  <- tree
          
          rules[[i]] <- rule
          # first base classifier's weight should be 1
          alphas[[i]] <- 1
          
          terms       <- tree$terms
          
          break
          
        }
        
        break
      }
      
      
      
      # Remove formulas since they waste memory
      if(i == 1){
        
        terms       <- tree$terms
        
      } else{
        
        tree$terms <- NULL
        
      }
      
      trees[[i]]  <- tree
      rules[[i]] <- rule
      
      alphas[[i]] <- alpha
      
      
    }
  }
  
  else{
    for(i in seq(n_rounds)){
      
      tree <- rpart::rpart(y ~ ., 
                           data = data.frame(X), weights = w,
                           method = "class", control = control, 
                           x = FALSE, y = FALSE, model = TRUE)
      
      rule <- rpart.rules(tree)
      pred <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
      
      # weighted error rate
      e    <- sum(w * (pred != y))
      
      if(e >= 0.5){
        
        e <- 1 - e
      }
      
      else{
        
        e <- e
        
      }
      # learning rate(weight) of each classifiers
      alpha <- 1/2 * log((1-e)/e)
      
      # update weight of each data
      w     <- w * exp(-alpha*pred*y)
      
      # normalize weight of each data
      w     <- w / sum(w)
      
      # If classifier's error rate is nearly 0, boosting process ends
      if(abs(e) < 1e-5){
        
        # if first base classifier predicts data perfectly, boosting process ends
        if(i == 1){
          
          trees[[i]]  <- tree
          
          rules[[i]] <- rule
          # first base classifier's weight should be 1
          alphas[[i]] <- 1
          
          terms       <- tree$terms
          
          break
          
        }
        
        break
      }
      
      
      
      # Remove formulas since they waste memory
      if(i == 1){
        
        terms       <- tree$terms
        
      } else{
        
        tree$terms <- NULL
        
      }
      
      trees[[i]]  <- tree
      rules[[i]] <- rule
      
      alphas[[i]] <- alpha
      
      
    }
  }
  
  
  
  
  out        <- list(alphas = unlist(alphas), 
                     trees = trees, 
                     terms = terms,
                     rules = rules)
  
  class(out) <- "adaboost"
  
  # create confusion matrix for in-sample fits
  y_hat                <- stats::predict(out, X)
  
  out$confusion_matrix <- table(y, y_hat)
  
  out
  
}

predict.adaboost <- function(object, X, type = c("response", "prob"),
                             n_tree = NULL, ...){
  # handle args
  type <- match.arg(type)
  
  if(is.null(n_tree)){
    
    tree_seq <- seq_along(object$alphas)
    
  } else{
    
    if(n_tree > length(object$alpha))
      
      stop('n_tree must be less than the number of trees used in fit')
    
    tree_seq <- seq(1, n_tree)
    
  }
  
  # evaluate score function on sample
  f <- 0
  
  for(i in tree_seq){
    
    tree       <- object$trees[[i]]
    
    tree$terms <- object$terms
    
    pred       <- as.integer(as.character(stats::predict(tree, data.frame(X),
                                                         type = "class")))
    f          <- f + object$alphas[i] * pred
    
  }
  
  # handle response type
  if(type == "response"){
    
    sign(f)
    
  } else if(type == "prob"){
    
    1/(1+exp(-2*f))
    
  }
  
}

kf <- function(tt){
  data    <- data
  
  # cross validation using adaboost to predict class of iris
  k       <- 5
  ##########ggogogogogogo
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  
  list    <- 1:k
  # data frame reset
  
  
  prediction_cm <- testset_copy_cm <- data.frame()
  prediction_bm <- testset_copy_bm <- data.frame()
  
  #function for k fold
  for(i in 1:k){
    
    # remove rows with id i from dataframe to create training set
    # select rows with id i to create test set
    trainset     <- subset(data, id %in% list[-i])
    
    testset      <- subset(data, id %in% c(i))
    
    #run a adaboost model
    model_cm        <- adaboost(trainset[, X_where], trainset[, y_where], TRUE , n_rounds = tt)
    model_bm      <- adaboost(trainset[, X_where], trainset[, y_where], FALSE , n_rounds = tt)  
    
    temp_cm         <- as.data.frame(predict(model_cm, testset))
    temp_bm         <- as.data.frame(predict(model_bm, testset))
    
    # 예측값을 예측 데이터 프레임의 끝에 추가.
    prediction_cm   <- rbind(prediction_cm, temp_cm)
    prediction_bm   <- rbind(prediction_bm, temp_bm)
    
    # 실제값(testset) testsetCopy에 추가.
    
    testset_copy_cm <- rbind(testset_copy_cm, as.data.frame(testset[, y_where]))
    testset_copy_bm <- rbind(testset_copy_bm, as.data.frame(testset[, y_where]))
    
    # 예측값과 실제값 데이터프레임.
    # add predictions and actual Sepal Length values
    result_cm            <- cbind(prediction_cm, testset_copy_cm[, 1])
    result_bm            <- cbind(prediction_bm, testset_copy_bm[, 1])
    
    
    names(result_cm)     <- c("Actual", "Predicted")
    names(result_bm)     <- c("Actual", "Predicted")
    
    confusion_matrix_cm  <- table(result_cm$Actual, result_cm$Predicted )
    confusion_matrix_bm  <- table(result_bm$Actual, result_bm$Predicted )
    
    accuracy_cm          <- sum(diag(confusion_matrix_cm)) / sum(confusion_matrix_cm)
    accuracy_bm          <- sum(diag(confusion_matrix_bm)) / sum(confusion_matrix_bm)
    
    result_cm            <- list("confusion_matrix_cm " = confusion_matrix_cm, "accuracy_cm" = accuracy_cm)
    
    result_bm            <- list("confusion_matrix_bm " = confusion_matrix_bm, "accuracy_bm" = accuracy_bm)
    
    
  }  
  temp1 <<- (result_cm$accuracy_cm)
  temp2 <<- (result_bm$accuracy_bm)
  cat(red(temp1))
  cat(("\n"))
  cat(blue(temp2))
  cat(("\n"))
}

gs <- function(qwe){
  for (i in 1:qwe){
    kf(i)
    cm_res <<- append(cm_res, temp1, after = length(cm_res))
    bm_res <<- append(bm_res, temp2, after = length(bm_res))
  }
  cm_res <<- as.vector(cm_res, mode = "numeric")
  bm_res <<- as.vector(bm_res, mode = "numeric")
  show()
}

show <- function(){
  cat(green("#############################\n"))
  cat(red(max(cm_res)))
  cat(("\n"))
  cat(blue(max(bm_res)))
  cat(("\n"))
  cat(red(mean(cm_res)))
  cat(("\n"))
  cat(blue(mean(bm_res)))
  plot(xx,cm_res, col = "red", ylim = c(min(min(bm_res), min(cm_res)),max(max(bm_res), max(cm_res))))
  lines(xx,cm_res,col = "red")
  par(new=TRUE)
  plot(xx,bm_res,col = "blue", xlab = "", ylab = "", axes = FALSE)
  lines(xx,bm_res,col = "blue")
}
####################################################################################################
####################################################################################################
gs(how_many)
