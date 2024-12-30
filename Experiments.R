# JR score

# Import the cleaned and imputed data 
data_THA <- read.csv("Data Import")


# Calculate the difference between pre and post-op HOOS JR score
Difference.HOOSJRScore <- data_THA$POST.HOOSJRScore - data_THA$PRE.HOOSJRScore

anch_MCID <- 17.7
dist_MCID <- 1/2 * sd(Difference.HOOSJRScore)
SCB <- 22

# Remove samples that cannot achieve any threshold due to a high pre-op score
tmp <- data_THA$PRE.HOOSJRScore < (100 - max(anch_MCID, dist_MCID, SCB )) # SCB
data_THA_2 <- data_THA[tmp, ]


#### matrix #####
### Univariate analysis with cat ####

X <- data_THA_2[,-which(colnames(data_THA_2) %in% c("MCID"))]
y <- data_THA_2[,"MCID"]
categ <- apply(X, 2, function(x){ # categorical variable = 2; binary variable = 1
  (length(names(table(x)))==2) + (length(names(table(x))) <5) ## binary: 2 ;categorical : 1
})
names(categ[categ==2])
names(categ[categ==1])
names(categ[categ==0])
X.bin <- X[,names(categ[categ==2])]
X.cat <- X[,names(categ[categ==1])]
X.num <- X[,names(categ[categ==0])]

# Binary variables
bin.chi <- lapply(colnames(X.bin), function(col.name){
  cont <- table(X.bin[,col.name], y)
  perc <- round(apply(cont, 2, function(x) x/sum(x) *100), 1)
  missing <- tapply(X.bin[,col.name], y, function(x) sum(is.na(x)))
  
  c(cont[2,], perc[2,], missing, round(chisq.test(cont)$p.value,3))
})


colnames(X.bin)
bin.chi.df <- do.call("rbind", bin.chi)
rownames(bin.chi.df) <- colnames(X.bin)
bin.chi.df
bin.chi.df2 <- data.frame(MCID0 = sprintf("%d (%.1f)", bin.chi.df[,1], bin.chi.df[,3]),
                          Missing0 = bin.chi.df[,5],
                          MCID1 =sprintf("%d (%.1f)", bin.chi.df[,2], bin.chi.df[,4]),
                          Missing1 = bin.chi.df[,6],
                          p.value = bin.chi.df[,7]
)
bin.chi.df
bin.chi.df2
rownames(bin.chi.df2)[bin.chi.df2$p.value<=0.10]

# Categorical variables
colnames(X.cat)
cat.chi <- lapply(colnames(X.cat), function(col.name){
  cont <- table(X.cat[,col.name], y)
  perc <- round(apply(cont, 2, function(x) x/sum(x) *100), 1)
  missing <- tapply(X.cat[,col.name], y, function(x) sum(is.na(x)))
  
  list(cbind(cont,perc), missing, round(chisq.test(cont)$p.value,3))
})

cat.chi.df <- do.call("rbind", lapply(cat.chi, `[[`, 1))
(chi.result <- cbind(do.call("rbind", lapply(cat.chi, `[[`, 2)),do.call("c", lapply(cat.chi, `[[`, 3))))
chi.result
which(chi.result[,3] <= 0.10)

cat.chi.df2 <- data.frame(MCID0 = sprintf("%d (%.1f)", cat.chi.df[,1], cat.chi.df[,3]),
                          MCID1 = sprintf("%d (%.1f)", cat.chi.df[,2], cat.chi.df[,4]))

cat.chi.df2
names(categ[categ==1])

colnames(X.cat)
row_names <- rep(colnames(X.cat), sapply(1:ncol(X.cat), function(idx){
  length(names(table(X.cat[,idx])))
}))
rownames(cat.chi.df2) <-paste(row_names,  rownames(cat.chi.df) ,sep="_")

# Continuous variables
con.chi <- lapply(colnames(X.num), function(col.name){
  dist <- tapply(X.num[,col.name], y, function(x) c(mean(x, na.rm=T), quantile(x, c(0.5, 0.25, 0.75), na.rm=T), sum(is.na(x))))
  test <- t.test(X.num[y==0, col.name], X.num[y==1, col.name], var.equal = ifelse(var.test(X.num[y==0, col.name], X.num[y==1, col.name])$p.value >0.05, T, F))
  c(do.call("c", dist), round(test$p.value,3))
  
})

con.chi.df <- do.call("rbind",con.chi)
rownames(con.chi.df) <- colnames(X.num)
con.chi.df
con.chi.df2 <- data.frame(MCID0 = sprintf("%.1f (%.1f, %.1f)", con.chi.df[,1], con.chi.df[,3], con.chi.df[,4]),
                          Missing0 = con.chi.df[,5],
                          MCID1 = sprintf("%.1f (%.1f, %.1f)", con.chi.df[,6], con.chi.df[,8], con.chi.df[,9]),
                          Missing1 = con.chi.df[,10],
                          p.value = con.chi.df[,11])

bin_con <- rbind(bin.chi.df2, con.chi.df2)

# Select significant variables
candidate_columns <- c(rownames(bin_con[bin_con$p.value<=.1,]), colnames(X.cat)[which(chi.result[,3] <= 0.10)])
candidate_columns


# Make a model matrix 
data_THA_2$thres <- anch_MCID # choose one of the threshold outcomes including anch_MCID, dist_MCID, and SCB
data_THA_2_mat <- model.matrix(~., data_THA_2)[,-1]

########### Experiments #############

## Significant variables - stepwise ###
mdl.glm <- glm(thres ~. ,data_THA_2, family="binomial")
summary(mdl.glm )

mdl.step <- step(mdl.glm)
(tmp <- summary(mdl.step))
exp(coef(mdl.step))


#### each label ####


### Penalized Logistic Regression ###

library(glmnet)
glmnet.mdl <- cv.glmnet(x = data_THA_2_mat[,-which(colnames(data_THA_2_mat) %in% c("thres"))],
                        y = data_THA_2_mat[, "thres"] , family = "binomial")
coef(glmnet.mdl, s = "lambda.min")

metrics <- function(test_label, test_vector){
  cont <- table(test_label, test_vector)
  recall <- ifelse(ncol(cont)>1, cont[2,2]/sum(cont[2,]), 0)
  prec <- ifelse(ncol(cont)>1, cont[2,2]/sum(cont[,2]), 0)
  return(c(
    acc = sum(diag(cont))/length(test_label),
    recall = recall,
    prec = prec,
    f1 = 2/(1/recall + 1/prec)
  ))
}

library(doParallel)
n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
repeated_glmnet <- foreach(repeated = 1:100, .packages = c("glmnet", "pROC", "PRROC"),.combine=rbind, .errorhandling = "remove") %dopar%{
  set.seed(repeated)
  test.idx <- sample.int(nrow(data_THA_2_mat), 0.3*nrow(data_THA_2_mat))
  train <- data_THA_2_mat[-test.idx, ]; test <- data_THA_2_mat[test.idx,]
  
  glmnet.mdl <- cv.glmnet(x = train[,-which(colnames(train) %in% c("MCID"))],
                          y = train[, "MCID"] , family = "binomial")
  
  scores_train <- predict(glmnet.mdl, newx = train[,-which(colnames(train) %in% c("MCID"))], s = "lambda.min", type ='response')
  
  perf <- do.call("rbind", lapply(seq(0,1,0.01), function(thres){
    c(thres = thres, metrics(train[,"MCID"], scores_train > thres))
  }))
  opt_thres <- seq(0,1,0.01)[which.max(perf[,"f1"])]
  
  scores_test <- predict(glmnet.mdl, newx= test[, -which(colnames(test) %in% c("MCID"))], s = "lambda.min", type = 'response')
  c(auc = auc(roc(test[,"MCID"], scores_test)), pr.integral = pr.curve(test[,"MCID"], scores_test)$auc.integral, 
    pr.dg = pr.curve(test[,"MCID"], scores_test)$auc.davis.goadrich, metrics(test[,"MCID"], scores_test > opt_thres))
  
}
stopCluster(cl)

res_glmnet <- c(colMeans(repeated_glmnet ),apply(repeated_glmnet, MARGIN = 2, function(x) quantile(x, c(.25,.75))))




### SVM ####
library(e1071); library(pROC); library(PRROC)

set.seed(1832)
test.idx <- sample.int(nrow(data_THA_2), 0.3*nrow(data_THA_2))
train <- data_THA_2[-test.idx, ]; test <- data_THA_2[test.idx,]

table(train$MCID);table(test$MCID)

svm.model <-  svm(as.factor(MCID) ~ ., data = train, probability=T)
summary(svm.model)
scores <- predict(svm.model, newdata = test[,-1],probability = T)
scores <- attr(scores, "probabilities")[,2]



repeated_svm_poly <- lapply(1:100, function(repeated){
  set.seed(repeated)
  test.idx <- sample.int(nrow(data_THA_2), 0.3*nrow(data_THA_2))
  train <- data_THA_2[-test.idx, ]; test <- data_THA_2[test.idx,]
  
  
  svm.model <-  svm(as.factor(MCID) ~ ., data = train, probability=T, kernel='polynomial', degree=3)
  scores_train <- predict(svm.model, newdata = train[,-1],probability = T)
  scores_train <- attr(scores_train, "probabilities")[,2]
  
  
  perf <- do.call("rbind", lapply(seq(0,1,0.01), function(thres){
    c(thres = thres, metrics(train$MCID, scores_train > thres))
  }))
  opt_thres <- seq(0,1,0.01)[which.max(perf[,"f1"])]
  
  scores_test <- predict(svm.model, newdata = test[,-1],probability = T)
  scores_test <- attr(scores_test, "probabilities")[,2]
  
  
  c(auc = auc(roc(test$MCID, scores_test)), pr.integral = pr.curve(test$MCID,scores_test)$auc.integral, 
    pr.dg = pr.curve(test$MCID, scores_test)$auc.davis.goadrich, metrics(test$MCID, scores_test > opt_thres))
})

repeated_svm_radial <- lapply(1:100, function(repeated){
  set.seed(repeated)
  test.idx <- sample.int(nrow(data_THA_2), 0.3*nrow(data_THA_2))
  train <- data_THA_2[-test.idx, ]; test <- data_THA_2[test.idx,]
  
  
  svm.model <-  svm(as.factor(MCID) ~ ., data = train, probability=T, kernel='radial', degree=3)
  scores_train <- predict(svm.model, newdata = train[,-1],probability = T)
  scores_train <- attr(scores_train, "probabilities")[,2]
  
  
  perf <- do.call("rbind", lapply(seq(0,1,0.01), function(thres){
    c(thres = thres, metrics(train$MCID, scores_train > thres))
  }))
  opt_thres <- seq(0,1,0.01)[which.max(perf[,"f1"])]
  
  scores_test <- predict(svm.model, newdata = test[,-1],probability = T)
  scores_test <- attr(scores_test, "probabilities")[,2]
  
  
  c(auc = auc(roc(test$MCID, scores_test)), pr.integral = pr.curve(test$MCID,scores_test)$auc.integral, 
    pr.dg = pr.curve(test$MCID, scores_test)$auc.davis.goadrich, metrics(test$MCID, scores_test > opt_thres))
})
repeated_svm_poly_df <- do.call('rbind',repeated_svm_poly)
res_svm_poly <- c(colMeans(repeated_svm_poly_df), apply(repeated_svm_poly_df,2, function(x) quantile(x, c(.25, .75))))
repeated_svm_radial_df <- do.call('rbind',repeated_svm_radial)
res_svm_radial <- c(colMeans(repeated_svm_radial_df), apply(repeated_svm_radial_df, 2, function(x) quantile(x, c(.25, .75))))


### RANDOM FOREST ####


generate_fold_number <- function(data, fold_no = 4){
  fold.number <- rep(1:fold_no, each = ceiling(nrow(data)/fold_no))
  fold.idx <- sample(fold.number, nrow(data))
  fold.idx
}


library(randomForest)
cv.rf <- function(data, parameter, fold_no = 5){
  sapply(parameter, function(ntree){
    fold.idx <- generate_fold_number(data, fold_no = fold_no)
    c(ntree=mean(sapply(1:fold_no, function(fn){
      train <- data[fold.idx != fn, ]
      val <- data[fold.idx == fn, ]
      
      mod_rf <- randomForest(x = train[,-which(colnames(train) %in% c("MCID"))], 
                             y = factor(train[,"MCID"]), mtry = mtry, ntree = ntree)
      
      p1 <- predict(mod_rf, val[,-which(colnames(val) %in% c("MCID"))], type = "prob")[,2]
      auc(roc(val$MCID, p1))
    })))
  })
}



parameters <- c(50, 100,150,200,250,300,350, 400)

library(doParallel)
n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(n_cores)

mtry <- 5; # options in this paper: 3, 4, and 5
registerDoParallel(cl)
repeated_cv_rf_5 <- foreach(repeated = 1:100, .packages = c("randomForest", "pROC", "PRROC"),.combine=rbind, .errorhandling = "remove") %dopar%{ # repeated_cv_rf_3, repeated_cv_rf_4, and repeated_cv_rf_5
  set.seed(repeated)
  test.idx <- sample.int(nrow(data_THA_2), 0.3*nrow(data_THA_2))
  train <- data_THA_2[-test.idx, ]; test <- data_THA_2[test.idx,]
  
  opt.par <- parameters[which.max(cv.rf(train, parameters, fold_no = 4))]
  opt.par
  
  
  mod_rf <- randomForest(x = train[,-which(colnames(train) %in% c("MCID"))], 
                         y = factor(train[,"MCID"]), mtry = mtry, ntree = opt.par)
  p_train <- predict(mod_rf, train[,-which(colnames(train) %in% c("MCID"))], type = "prob")[,2]
  
  perf <- do.call("rbind", lapply(seq(0,1,0.01), function(thres){
    c(thres = thres, metrics(train$MCID, p_train > thres))
  }))
  opt_thres <- seq(0,1,0.01)[which.max(perf[,"f1"])]
  
  p_test <- predict(mod_rf, test[,-which(colnames(test) %in% c("MCID"))], type = "prob")[,2]
  c(auc = auc(roc(test$MCID, p_test)), pr.integral = pr.curve(test$MCID, p_test)$auc.integral, 
    pr.dg = pr.curve(test$MCID, p_test)$auc.davis.goadrich, metrics(test$MCID, p_test > opt_thres))
}
stopCluster(cl)

res_rf_3 <- c(colMeans(repeated_cv_rf_3), apply(repeated_cv_rf_3, 2, function(x) quantile(x, c(.25, .75))))
res_rf_4 <- c(colMeans(repeated_cv_rf_4), apply(repeated_cv_rf_4, 2, function(x) quantile(x, c(.25, .75))))
res_rf_5 <- c(colMeans(repeated_cv_rf_5), apply(repeated_cv_rf_5, 2, function(x) quantile(x, c(.25, .75))))


rbind(res_glmnet,
      res_svm_poly ,res_svm_radial,
      res_rf_3 ,res_rf_4 ,res_rf_5) 