

p <- 1000

X <- c( rnorm(100, -2, 1), rnorm(100, 2, 1), rnorm(800, 0, ) )

q <- 0.05

## Right-sided test
p.right <- pnorm( X, lower.tail=FALSE )
clat.right <- CLAT( p.right, q)
rejind <- which( clat.right$SigInd!=0 )

## Right-sided test
p.left <- pnorm( X, lower.tail=TRUE )
clat.left <- CLAT( p.left, q)
rejind <- which( clat.left$SigInd!=0 )


## Two-sided test
p.right <- pnorm( X, lower.tail=FALSE )
clat.right <- CLAT( p.right, q)

p.left <- pnorm( X, lower.tail=TRUE )
clat.left <- CLAT( p.left, q)
rejind <- c( which( clat.left$SigInd!=0 ), which(clat.right$SigInd!=0) ) 
