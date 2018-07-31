### MNIST EXAMPLE ###

to.read = file("C:/users/gravesee/Downloads/t10k-images.idx3-ubyte", "rb")
readBin(to.read, integer(), n=4, endian="big")

images <- sapply(seq.int(10000), function(x) {
  readBin(to.read,integer(), size=1, n=28*28, endian="big")
})

par(mfrow=c(5, 4))
par(mar=c(0,0,0,0))
for (i in 1:20) image(matrix(images[,i] < 0, 28, 28)[,28:1], axes=F)
par(mfrow=c(1, 1))

library(Matrix)
d <- t(images)
d <- Matrix((d < 0) + 0L)
z <- as(as(d, "nsparseMatrix"), "ngCMatrix")




set.seed(100)
res <- BMM(z, K=20L, max.iter = 20L, verbose = 1L, hbbmm = 1L)

res <- BMM(z, K=20, max.iter = 20, verbose = TRUE, hbbmm = TRUE)

p <- predict(res, z[1:100,])

z <- res$prototypes

par(mfrow=c(5,4))
par(mar=c(0,0,0,0))
for (i in seq.int(nrow(z))) {
  image(matrix(res$prototypes[i,], 28, 28)[,28:1], axes=F)
}
par(mfrow=c(1,1))
