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

res <- BMM(z, K=10L, max.iter = 25L, verbose = 1L)

z2 <- matrix(as.integer(as.matrix(d)), nrow(d), ncol(d))


p <- predict_dense_matrix(z2, res$prototypes, res$pis)

p <- predict_sparse_matrix(z, res$prototypes, res$pis)



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

iso <- iForest(titanic[-1], nt = 100, phi=8)
nodes <- predict(iso, titanic[-1], sparse=T)
nodes <- as(as(nodes, "nsparseMatrix"), "ngCMatrix")

bmm <- BMM(nodes, K=10L, max.iter = 100L, verbose = 1L)

p <- predict_sparse_matrix(nodes, bmm$prototypes, bmm$pis)

library(onyx)

mod <- bin(data.frame(p$z), y=titanic$Survived)

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

