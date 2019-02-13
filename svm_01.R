library("quadprog")
library("Matrix")
library("Distance")
library("rdist")
library("lattice") # heatmap
library("MASS") # heatmap color



iris <- read.csv(file = "iris_data.csv", header = TRUE)


n_class = (0:2)    # class ÀÇ °³¼ö 0,1,2
var_data = (1:4)   # variable

data_train_num <- sample(1:nrow(iris), size = round(0.6 * nrow(iris))) # train ratio
data_train = iris[data_train_num,] # train data
data_test = iris[-data_train_num,] # test data

svm_cal <- function(p1,p2){ # par_1, par2_
  res_svm <- list()
  res = 0
  for (ii in n_class){
    data_train$newclass <- ifelse(data_train$class == ii, 1, -1)
    data_test$newclass <- ifelse(data_test$class == ii, 1, -1)
    
    euc_matrix_1 <- t(t(data_train$newclass)) %*% t(data_train$newclass)                                                           # y y^T
    euc_matrix_2 <- exp((as.matrix(dist(data_train[var_data], method = "euclidean", diag = T, upper = T)**2))*(-p1))            # k(xi,xj)
    
    
    Dmat <- as.matrix((nearPD(euc_matrix_1 * euc_matrix_2))$mat)                                                              # Dmat
    
    dvec <- as.matrix(rep(1, nrow(data_train)))                                                                                          # dvec
    
    Amat <- t(rbind(t(data_train$newclass), rbind(diag(1, nrow(data_train), nrow(data_train))), diag(-1, nrow(data_train), nrow(data_train))))
    bvec <- rbind((t(t(rep(0,1)))), (rbind(t(t(rep(0, nrow(data_train)))), t(t(rep(-p2, nrow(data_train)))))))
    
    
    alphas <- solve.QP(Dmat, dvec, Amat, bvec)$solution
    
    bias <- mean(t(t(data_train$newclass)) - (diag(c(t((diag(alphas) %*% t(t(data_train$newclass)))))) %*% (Dmat %*% t(t(dvec)))))
    xxx <- list()
    
    for (k1 in 1:nrow(data_test)){
      for (k2 in 1:nrow(data_train)){
        ttt <- exp((cdist(data_test[var_data][k1,], data_train[var_data][k2,], metric = "euclidean", p = 2)**2)*(-p1))
        xxx <- append(xxx, ttt, after = length(xxx))
      }
    }
    
    m2 <- matrix(unlist(xxx), nrow = nrow(data_test), ncol = nrow(data_train), byrow = TRUE)
    m3 <-(alphas * t(data_train$newclass)) #a1y1 a2y2 
    m4 <- do.call(rbind, replicate(nrow(data_test), m3, simplify=FALSE)) # a1y1 a2y2 many
    m5 <- m4 * m2
    m6 <- m5 %*% t(t(dvec))
    
    tt = 0
    
    for (i in 1:nrow(m6)){
      if (m6[i] > 0 & data_test$newclass[i] == 1){
        tt = tt + 1
      }
      else if (m6[i] < 0 & data_test$newclass[i] == -1){
        tt = tt + 1
      }
    }
    
    res = res + tt * 100 / nrow(data_test)
  }
return (res / length(n_class))
}
#svm_cal(0.001, 0.001)

svm <- function(m1,m2,m3,m4){ # min max boundary of (rbf gamma / c) input 
  c1 <- log(m2) - log(m1)
  c2 <- log(m4) - log(m3)
  mat_res <- matrix(nrow = c1, ncol = c2)
  list_res <- list()
  list_gamma <- list()
  list_c <- list()
  
  for (i in 1:c1){
    list_gamma <- append(list_gamma, m1 * 10^(i - 1), after = length(list_gamma))
  }
  
  for (i in 1:c2){
    list_c <- append(list_c, m1 * 10^(i - 1), after = length(list_c))
  }
  
  #mat_gamma <-do.call(rbind, replicate(length(list_c), as.matrix(unlist(list_gamma)), simplify = FALSE))
  #mat_C <-do.call(rbind, replicate(length(list_gamma), as.matrix(unlist(list_c)), simplify = FALSE))
  
  for (i1 in 1:c1){
    for (i2 in 1:c2){
      list_res <- append(list_res, svm_cal(m1 * 10^(i1 - 1), m3 * 10^(i2 - 1)), after = length(list_res))
      }
    }
  z_mat <- matrix(unlist(list_res), nrow = length(list_gamma), ncol = length(list_c), dimnames = list(c(list_gamma), c(list_c)), byrow = TRUE)
  
  print(z_mat)
  #heatmap(z_mat, Rowv = rownames(z_mat), Colv = rownames(z_mat), col = cm.colors(256), scale = "column", margins = c(5,10), dendrogram = "none")
  #heatmap.2(z_mat, Rowv = rownames(z_mat), Colv = colnames(z_mat), dendrogram = 'none')
  #matrix.heatmap(z_mat, Rowv(rownames(z_mat)), Colv(colnames(z_mat)))
  levelplot(z_mat)
  print(max(z_mat))
  print (which(z_mat == max(z_mat), arr.ind = TRUE))
  
  }





svm(0.01,100,0.01,100) # 0.001 under error??







