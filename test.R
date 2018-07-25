to.read = file("C:/users/gravesee/Downloads/t10k-images.idx3-ubyte", "rb")
readBin(to.read, integer(), n=4, endian="big")

images <- sapply(seq.int(10000), function(x) {
  readBin(to.read,integer(), size=1, n=28*28, endian="big")
})

l <- file("c:/Users/gravesee/Downloads/t10k-labels.idx1-ubyte", "rb")
readBin(l, integer(), size = 4, n=2, endian = "big")
labels <- readBin(l, integer(), size = 1, n=10000, endian = "big")

d <- t(images)
d <- Matrix((d < 0) + 0L)
z <- as(as(d, "nsparseMatrix"), "ngCMatrix")

set.seed(100)
res1 <- BMM(z, K=10L, max.iter = 100L, verbose = 1L, hbbmm = 0L)

set.seed(100)
res2 <- BMM(z, K=10L, max.iter = 100L, verbose = 1L, hbbmm = 1L)


z2 <- matrix(as.integer(as.matrix(d)), nrow(d), ncol(d))


p <- predict_dense_matrix(z2, res$prototypes, res$pis)

p <- predict_sparse_matrix(z, res$prototypes, res$pis)


z <- res$prototypes

par(mfrow=c(4,3))
par(mar=c(0,0,0,0))
for (i in seq.int(10)) {
  image(matrix(res2$prototypes[i,], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))

## example images
par(mfrow=c(5,4))
par(mar=c(0,0,0,0))
for (i in seq.int(20)) {
  image(matrix(d[i,], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))




fit <- MASS::lda(p$z, grouping=labels)

preds <- predict(fit)

hist(p$z[,5])






### plot images

par(mfrow=c(4,3))
par(mar=c(0,0,0,0))
for (i in seq.int(10)) {
  image(matrix(res$prototypes[i,], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))

sapply(split(labels, res$cluster), table)

cl17 <- which(res$cluster == 17)
z <- t(images)
image(matrix(z[cl17[3],], 28, 28)[,28:1], axes=F)

par(mfrow=c(5,5))
par(mar=c(0,0,0,0))
for (i in seq.int(25)) {
  image(matrix(z[i,], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))



data(titanic, package="onyx")
library(isofor)

iso <- iForest(titanic[-1], nt = 20, phi=4)
nodes <- predict(iso, titanic[-1], sparse=T)
nodes <- as(as(nodes, "nsparseMatrix"), "ngCMatrix")

f <- titanic$Survived == 1
bmm_1 <- BMM(nodes[f,], K=1L, max.iter = 100L, verbose = 1L)
bmm_0 <- BMM(nodes[!f,], K=1L, max.iter = 100L, verbose = 1L)

## combine the prototypes from each

prototypes <- rbind(bmm_1$prototypes, bmm_0$prototypes)
pis <- rep(1, nrow(prototypes))/nrow(prototypes)

z <- predict_sparse_matrix(nodes, prototypes, pis)

head(z$z)

pred <- apply(z$z, 1, which.max)

table(titanic$Survived, pred %in% 1:3)

library(ks)

## create pairwise scores
scores <- data.frame(score=-z$z[,1] + z$z[,2], perf=titanic$Survived)
scores <- cbind(scores, z$z)
names(scores)[3:4] <- c("z1", "z2")

tbl <- ks_table(perf~score+z1+z2, data=scores, number_bins = 10)

scores$perf <- titanic$Survived
scores$diff_14 <- scores$X1 - scores$X4
scores$diff_15 <- scores$X1 - scores$X5
scores$diff_16 <- scores$X1 - scores$X6

scores$avg_1s <- rowMeans(z$z[,1:3])
scores$avg_0s <- rowMeans(z$z[,-(1:3)])

diffs <- list()
idx <- 1
for (i in 1:3) {
  for (j in 4:8) {
    diffs[[idx]] <- z$z[,i] - z$z[,j]
    names(diffs)[idx] <- paste0("diff_", paste0(i, j))
    idx <- idx + 1
  }
}

scores <- data.frame(diffs)
scores$perf <- titanic$Survived
scores <- cbind(scores, z$z)

tbls <- ks_table(perf~diff_38, data=scores, number_bins = 10)

library(onyx)
mod <- bin(scores, scores$perf, min.cnt=25, min.res=10, mono=2)

mod <- bin(data.frame(z$z), scores$perf, min.cnt=25, min.res=10, mono=2)

mod$sort()

score <- z$z[,3] - z$z[,7]

library(ggplot2)

plt <- data.frame(score, Survived=titanic$Survived)
ggplot(plt, aes(x=score, y=Survived)) + geom_smooth(method="glm", method.args=list(family="binomial"))

library(pROC)

r1 <- pROC::roc(plt$Survived, plt$score)
plot(r1)


## predict two BMMs, one on survived one not
bmm1






## diffs
l <- list()
idx <- 1
for (i in seq.int(ncol(z$z) - 1)) {
  for (j in (i+1):ncol(z$z)) {
    l[[idx]] <- z$z[,i] - z$z[,j]
    idx <- idx + 1
  }
}

X <- data.frame(l)

z <- predict_sparse_matrix(nodes, bmm$prototypes, bmm$pis)



library(onyx)
names(X) <- make.names(seq.int(ncol(X)))
mod <- bin(X, y=titanic$Survived, mono=2, min.res = 10)
mod$sort()
mod$set_step(lvl=1)

tapply(titanic$Survived, bmm$cluster, mean)


library(spamming)

dst <- as.dist(spamming(nodes))
cl <- hclust(dst)



P1 <- c(0.9, 0.9, 0.9, 0.1, 0.1)
P2 <- c(0.1, 0.1, 0.9, 0.9, 0.9)

prototypes <- list(P1, P2)
weights <- c(0.25, 0.75)

x <- t(replicate(1000, {
  
  ## pick a random prototype
  i <- sample(1:2, size = 1, prob =  weights)
  
  ## sample bits from the chosen prototype
  sapply(prototypes[[i]], function(p) rbinom(1, 1, p))
  
}))


## Training a BMM 
set.seed(1234)
res <- BMM(data = x, K = 2L, max.iter = 20L, verbose = 1L)


predict_dense_matrix(x, t(res$prototypes), res$pis)

## logsumexp









## chop up the titanic dataset
f <- function(x) {
  if (is.numeric(x)) {
    qs <- quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE)
    res <- cut(x, c(-Inf, unique(qs), Inf))
    levels(res) <- c(levels(res), "MISSING")
    res[is.na(x)] <- "MISSING"
    res
  } else {
    x
  }
}

library(onyx)
data(titanic)

xs <- lapply(titanic[-1], f)
X <- data.frame(xs)
m <- Matrix::sparse.model.matrix(~., X)
M <- as(m, "nsparseMatrix")

library(isofor)
iso <- iForest(titanic[-1], nt=50, phi = 8)
nodes <- predict(iso, titanic[-1], sparse = TRUE)

M <- as(nodes, "nsparseMatrix")
M <- as(M, "ngCMatrix")

bmm <- BMM(M, K=10L, max.iter = 100L, verbose = 1L)

pro <- bmm$prototypes[c(1,4:5),]

z <- predict_sparse_matrix(M, pro, c(0.5, 0.5, 0.5))

score <- z$z[,3] - z$z[,1]
library(pROC)

r1 <- roc(titanic$Survived, score)
r1 <- roc(titanic$Survived, z$z[,3])

mod <- bin(data.frame(z$z), titanic$Survived, mono=2, min.res=10, min.cnt=25)






### Bayes factors for cluster 1





