# ===========================================
# |||||||||||||| LOAD DATA ||||||||||||||||||
# ===========================================
data<-foreign::read.arff("11.arff")
col_names<-colnames(data)
col_names<-make.names(col_names, uniq=TRUE)
colnames(data)<-col_names

# ===========================================
# ||||||||||||| NOT SUPERVISED ||||||||||||||
# ===========================================

library(mclust)
library(gridExtra)
library(factoextra)

# Pre-processing
data.all <- data
row.names <- paste(data.all$class, rownames(data.all), sep = "." ); row.names[1:4]
data.all$class <- NULL
data.all <- scale(data.all, center = TRUE, scale = TRUE)
rownames(data.all) <- row.names

# k-means
p1 <- fviz_nbclust(data.all, kmeans, method = "wss")
p2 <- fviz_nbclust(data.all, kmeans, method = "silhouette")

km.2 <- kmeans(data.all, 2, nstart = 10)
p3 <- fviz_cluster(km.2, data.all, ellipse.type = "norm",
                   ellipse.level = 0.95,
                   ellipse.alpha = 0,
                   main = "K-Means (k = 2)",
                   geom = c("point"), 
                   ggtheme =   theme_bw(base_family = "Times New Roman",
                                        base_size = 15)
)

km.4 <- kmeans(data.all, 4, nstart = 10)
p4 <- fviz_cluster(km.4, data.all, ellipse.type = "norm",
                   ellipse.level = 0.95,
                   ellipse.alpha = 0,
                   main = "K-Means (k = 4)",
                   geom = c("point"), 
                   ggtheme =   theme_bw(base_family = "Times New Roman",
                                        base_size = 15)
)


grid.arrange(p1, p2, p3, p4, nrow = 2)



# ===========================================
# ||||||||||||||| RELIEF FFS ||||||||||||||||
# ===========================================
library(mlr)
task <- makeClassifTask(id = "data", data = data, target = "class", positive = "Tumor")
task  # To see task information
norm.task <- normalizeFeatures(task, method = "standardize")  # (X - MEAN) / STD

get_top <-function(x, n){
  best <- x$data
  best <- best[order(-best$value), ]
  best <- best[1:n, ]
  return(best)
}

fv <- generateFilterValuesData(norm.task, method = c("FSelector_relief"))   # Multivariate (take time)
top.relief <- get_top(fv, 100)

top.relifef.name <- top.relief$name
top.data <- data[,top.relifef.name]
top.data$class <- data$class

#write.csv(top.data, file = "RELIEF_TOP_miRNA.csv")   # To save the results

RELIEF <- read.table("RELIEF_TOP_miRNA.csv", header = T, sep = ",", row.names = 1)


# Select RELIEF USING A SEQUENTIAL SEARCH + SVM


relief.task <- makeClassifTask(id = "RELIEF miRNA", 
                               data = RELIEF, target = "class",
                               positive = "Tumor")

relief.task <- normalizeFeatures(relief.task, method = "standardize")


sfeats.svm <- selectFeatures(learner = makeLearner("classif.svm", kernel = "radial"),
                             task = relief.task,
                             control = makeFeatSelControlSequential(maxit = 200,
                                                                    method = 'sfs',
                                                                    max.features = 15),
                             resampling = makeResampleDesc("LOO"), # LOO
                             measures = list(acc))
                            
sfeats.svm



sfeats.knn <- selectFeatures(learner = makeLearner("classif.kknn"),
                             task = relief.task,
                             control = makeFeatSelControlSequential(maxit = 200,
                                                                    method = 'sfs',
                                                                    max.features = 15),
                             resampling = makeResampleDesc("LOO"), # LOO
                             measures = list(acc))
sfeats.knn


most.important.miRNA<- unique(c(sfeats.svm$x, sfeats.knn$x))
most.important.miRNA <- c("miR.183.", "miR.99b", "miR.497.", "miR.139.3p","Candidate.12.3p")
# THIS miRNA HAS BEEN REMOVED -> Candidate.12.5p (0 COUNTS IN ALL INDIVIDUALS (ONE OR TWO HAS 1 COUNT, NOT SIGNIFICANT))
labels <- RELIEF$class
RELIEF <- RELIEF[,most.important.miRNA]
RELIEF$class <- labels


# ===========================================
# |||||||||||| SUPERVISED CLAS.||||||||||||||
# ===========================================
norm.task # all miRNA

relief.task <- makeClassifTask(id = "RELIEF miRNA", 
                               data = RELIEF, target = "class",
                               positive = "Tumor")

relief.task <- normalizeFeatures(relief.task, method = "standardize")

rdesc <- makeResampleDesc("RepCV", folds = 10, stratify = T)  # 10 Repeated 10-CV
 

# SVM ========================================

res.svm.all <- resample(learner = makeLearner("classif.svm", 
                                              kernel = "radial",
                                              id = "SVM.Original_miRNA",
                                              predict.type = "prob"),
                        task = norm.task,
                        resampling = rdesc,
                        measures = list(acc, tpr, tnr))

# Aggregated Result: acc.test.mean=0.7625000,tpr.test.mean=0.8500000,tnr.test.mean=0.6783333

res.svm.relief <- resample(learner = makeLearner("classif.svm", 
                                                 kernel = "radial",
                                                 id = "SVM.miRNA_Selection",
                                                 predict.type = "prob"),
                           task = relief.task,
                           resampling = rdesc,
                           measures = list(acc, tpr, tnr))

# Aggregated Result: acc.test.mean=0.8523333,tpr.test.mean=0.9116667,tnr.test.mean=0.7950000







# wKNN =========================================

res.knn.all <- resample(learner = makeLearner("classif.kknn", 
                                              id = "KNN.Original_miRNA",
                                              predict.type = "prob"),
                        task = norm.task,
                        resampling = rdesc,
                        measures = list(acc, tpr, tnr))
# Aggregated Result: acc.test.mean=0.7840000,tpr.test.mean=0.9000000,tnr.test.mean=0.6650000

res.knn.relief <- resample(learner = makeLearner("classif.kknn", 
                                                 id = "KNN.miRNA_Selection",
                                                 predict.type = "prob"),
                           task = relief.task,
                           resampling = rdesc,
                           measures = list(acc, tpr, tnr))
# Aggregated Result: acc.test.mean=0.9358333,tpr.test.mean=0.9650000,tnr.test.mean=0.9050000


df.p <- generateThreshVsPerfData(list(res.svm.all, res.svm.relief, res.knn.all, res.knn.relief), 
                                 measures = list(fpr, tpr, mmce))
plotROCCurves(df.p)

df.n <- generateThreshVsPerfData(list(res.svm.all, res.svm.relief, res.knn.all, res.knn.relief), 
                                 measures = list(fnr, tnr, mmce))
plotROCCurves(df.n)


# PERFORMANCE ANALYSIS ========================

# BOXPLOT ACCURACY ============================
svm.all.acc <- res.svm.all$measures.test$acc
svm.relief.acc <- res.svm.relief$measures.test$acc
knn.all.acc <- res.knn.all$measures.test$acc
knn.relief.acc <- res.knn.relief$measures.test$acc

accuracy <- c(svm.all.acc, svm.relief.acc, knn.all.acc, knn.relief.acc)
ALGORITHM <- c(rep("SVM", 200), rep("KNN", 200))
probe <- c(rep("SVM.ALL", 100),rep("SVM.RELIEF+SEQ", 100),
           rep("KNN.ALL", 100),rep("KNN.RELIEF+SEQ", 100))


performance <- data.frame(accuracy, ALGORITHM, probe)

plot <- ggplot(performance, aes(x = probe, y = accuracy, fill = ALGORITHM)) + 
  scale_fill_manual(values=c("#000000", "#D8D8D8"), name = "Classifier") +  # To select class color
  geom_boxplot(alpha = 0.5, colour = "black") + 
  scale_x_discrete(name = "Algorithm + Feature selection") +
  scale_y_continuous(name = "Accuracy 10 Repeated 10-CV", breaks = round(seq(0.0, 1, by = 0.1), 2),
                     limits = c(0.0,1)) +                                # 5 -CV
  ggtitle("Supervised Classification") +
  theme(plot.title = element_text(hjust = 0.6, size = 20, family = 'Times New Roman'),
        panel.border = element_blank(),
        panel.grid.minor = element_line(colour = "grey80"),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black")) 
plot 

# BARCHART MEAN FNR/FPR ============================
#FPR
svm.all.tpr <- mean(res.svm.all$measures.test$tpr); svm.all.tpr
svm.relief.tpr <- mean(res.svm.relief$measures.test$tpr); svm.relief.tpr
knn.all.tpr <- mean(res.knn.all$measures.test$tpr); knn.all.tpr
knn.relief.tpr <- mean(res.knn.relief$measures.test$tpr); knn.relief.tpr

label <- rep("TPR", 4)
tpr <- c(svm.all.tpr, svm.relief.tpr,knn.all.tpr, knn.relief.tpr )
true.positive.rate <- data.frame(tpr, label, c("SVM.ALL", "SVM.RELIEF+SEQ", "KNN.ALL", "KNN.RELIEF+SEQ"))
colnames(true.positive.rate) <- c("values", "label", "probe")

#FNR
svm.all.tnr <- mean(res.svm.all$measures.test$tnr);svm.all.tnr
svm.relief.tnr <- mean(res.svm.relief$measures.test$tnr); svm.relief.tnr
knn.all.tnr <- mean(res.knn.all$measures.test$tnr);knn.all.tnr
knn.relief.tnr <- mean(res.knn.relief$measures.test$tnr); knn.relief.tnr

label <- rep("TNR", 4)
tnr <- c(svm.all.tnr, svm.relief.tnr,knn.all.tnr, knn.relief.tnr )
true.negative.rate <- data.frame(tnr, label, c("SVM.ALL", "SVM.RELIEF+SEQ", "KNN.ALL", "KNN.RELIEF+SEQ"))
colnames(true.negative.rate) <- c("values", "label", "probe")



t.positive.negative <- rbind(true.positive.rate, true.negative.rate)


ggplot(t.positive.negative, aes(x = probe, y = values)) + 
  scale_fill_manual(values=c("#000000", "#666E67"), name = "TNR/TPR") +  # To select class color
  geom_bar(aes(fill = label), stat = "identity",position = "dodge", alpha = 0.8) +
  scale_x_discrete(name = "Mean true positive and negative rates") +
  scale_y_continuous(name = "True negative/positive rate", breaks = round(seq(0,1, by = 0.1), 2),
                     limits = c(0,1)) + 
  ggtitle("True positive & True negative rates") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, family = 'Times New Roman'),
        panel.border = element_blank(),
        panel.grid.minor = element_line(colour = "grey80"),
        panel.background = element_rect(fill = "white", colour = "white"),
        
        axis.line = element_line(colour = "black")) 





# ===========================================
# ||||||||||||||||||| PCA |||||||||||||||||||
# ===========================================
library(matlib)
library(maptools)
library(RColorBrewer)
library(ggfortify)
# PCA ANALYSIS ==============================
PCA <- princomp(~miR.183.+miR.99b+miR.497.+miR.139.3p+Candidate.12.3p, cor=TRUE, data=RELIEF)

# Loadings
print(unclass(loadings(PCA)))

# Variances
print(PCA$sd^2)

# Contribution of each feature to the components (loadings)
plot(PCA$loadings, xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), xlab="Principal Component 1", ylab="Principal Component 2",type="n") # display the graphic
text(PCA$loadings, labels=row.names(PCA$loadings),cex=1,pos=c(4,4,4,2,2,4))
PCA$loadings

# Selected features
autoplot(prcomp(RELIEF[,1:5]), data = RELIEF, colour = 'class')


# ===========================================
# ||||||||||||| NOT SUPERVISED ||||||||||||||
# ===========================================

library(mclust)
library(gridExtra)
library(factoextra)

# Pre-processing
data.all <- RELIEF
row.names <- paste(data.all$class, rownames(data.all), sep = "." ); row.names[1:4]
data.all$class <- NULL
data.all <- scale(data.all, center = TRUE, scale = TRUE)
rownames(data.all) <- row.names

# k-means
p1 <- fviz_nbclust(data.all, kmeans, method = "wss")
p2 <- fviz_nbclust(data.all, kmeans, method = "silhouette")

km.2 <- kmeans(data.all, 2, nstart = 20)
p3 <- fviz_cluster(km.2, data.all, ellipse.type = "norm",
                   ellipse.level = 0.95,
                   ellipse.alpha = 0,
                   main = "K-Means (k = 2)",
                   geom = c("point"), 
                   ggtheme =   theme_bw(base_family = "Times New Roman",
                                        base_size = 15)
)

km.3 <- kmeans(data.all, 4, nstart = 10)
p4 <- fviz_cluster(km.3, data.all, ellipse.type = "norm",
                   ellipse.level = 0.95,
                   ellipse.alpha = 0,
                   main = "K-Means (k = 4)",
                   geom = c("point"), 
                   ggtheme =   theme_bw(base_family = "Times New Roman",
                                        base_size = 15)
)


grid.arrange(p1, p2, p3, p4, nrow = 2)


cluster.1 = km.2$cluster[ km.2$cluster == 1]   # 8 Normal
cluster.2 = km.2$cluster[ km.2$cluster == 2]   # 21 Normal 29 Tumor
cluster.1
cluster.2

cluster.1= km.3$cluster[ km.3$cluster == 1]    
cluster.2= km.3$cluster[ km.3$cluster == 2]    
cluster.3= km.3$cluster[ km.3$cluster == 3]    
cluster.4= km.3$cluster[ km.3$cluster == 4]    
cluster.1   # 6 NORMAL
cluster.2   # 1 NORMAL
cluster.3   # 3 TUMOR
cluster.4   # 22 NORMAL / 26 TUMOR


#=============================================================
# ||||||||||||||||||| GAUSSIAN MIXTURES ||||||||||||||||||||||
#=============================================================
BIC <- mclustBIC(data.all)
plot(BIC)
summary(BIC)

mod1 <- Mclust(data.all, x = BIC)
plot(mod1, what = "classification")
table(labels, mod1$classification)

#labels    1  2  3  4  5  6  7
#Normal  4 17  6  1  0  0  1
#Tumor   1  7  0  0 12  3  6


#=============================================================
# ||||||||||||||||||| HISTOGRAMS ||||||||||||||||||||||
#=============================================================
library(ggplot2)
library(gridExtra)
library(grid)

df <- read.table("SELECTED_MIRNA.csv", header = T, row.names = 1, sep = ",")

miR.183. <- df[,c("miR.183.","CONDITION")]
miR.99b <- df[,c("miR.99b","CONDITION")]
miR.497. <- df[,c("miR.497.","CONDITION")]
miR.139.3p <- df[,c("miR.139.3p","CONDITION")]
Candidate.12.3p <- df[,c("Candidate.12.3p","CONDITION")]


p1 <- ggplot(miR.183., aes(x=miR.183., fill = CONDITION)) + stat_bin(bins = 30) +
  scale_fill_manual(values=c("#FF6347", "#00BFFF"), name = "") +
  geom_histogram(alpha = 0.2) +  theme_classic() + theme(legend.position="top") 
  

p2 <- ggplot(miR.99b, aes(x=miR.99b, fill = CONDITION)) + stat_bin(bins = 30) +
  geom_histogram(alpha = 0.5) +  scale_fill_manual(values=c("#FF6347", "#00BFFF"), name = "") +
  theme_classic() + theme(legend.position="top")

p3 <- ggplot(miR.497., aes(x=miR.497., fill = CONDITION)) + stat_bin(bins = 30) +
  geom_histogram(alpha = 0.5) +   scale_fill_manual(values=c("#FF6347", "#00BFFF"), name = "") +
  theme_classic() + theme(legend.position="top") 

p4 <- ggplot(miR.139.3p, aes(x=miR.139.3p, fill = CONDITION)) + stat_bin(bins = 30) +
  geom_histogram(alpha = 0.5) +  scale_fill_manual(values=c("#FF6347", "#00BFFF"), name = "") +
  theme_classic() + theme(legend.position="top")

p5 <- ggplot(Candidate.12.3p, aes(x=Candidate.12.3p, fill = CONDITION)) + stat_bin(bins = 30) +
  geom_histogram(alpha = 0.5) +   scale_fill_manual(values=c("#FF6347", "#00BFFF"), name = "") +
  theme_classic() + theme(legend.position="top")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3)


