to.read = file("~/Downloads/t10k-images-idx3-ubyte", "rb")
readBin(to.read, integer(), n=4, endian="big")

images <- sapply(seq.int(10000), function(x) {
  readBin(to.read,integer(), size=1, n=28*28, endian="big")
})

l <- file("~/Downloads/t10k-labels-idx1-ubyte", "rb")
readBin(l, integer(), size = 4, n=2, endian = "big")
labels <- readBin(l, integer(), size = 1, n=10000, endian = "big")

d <- t(images)
d <- (d > 0) + 0L

res <- BMM(d, K=16L, max.iter = 100L, verbose = 1L)


### plot images

par(mfrow=c(4,4))
par(mar=c(0,0,0,0))
for (i in seq.int(16)) {
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

