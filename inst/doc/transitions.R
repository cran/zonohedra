## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_opt = options( width=144 )

# if( !file.exists("figs") ) dir.create("figs")

require("rgl",quietly=TRUE)
rgl::setupKnitr(autoprint = TRUE)

## ----echo=TRUE,  message=FALSE----------------------------------------------------------------------------------------------------------------
library(zonohedra)

## ----echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6.5, fig.height=4, fig.cap='Figure 2.1  four points in the 2-transition complex, visualized with bar graphs', out.width="100%", cache=FALSE----
mybarplot <- function( x )  {
n = length(x)
plot( c(0,n), c(0,1), type='n', tcl=0, las=1, xaxt='n', xlab='', ylab='', mgp=c(3,0.25,0) )
grid( nx=NA, ny=NULL, lty=1 )
barplot( x, names.arg=1:n, space=0, add=T, yaxt='n', mgp=c(3,0.25,0) )
}

x1 = numeric(10) ; x1[ c(3,8) ] = exp( c(-0.25,-1) ) ; x1[ 4:7 ] = 1
x2 = numeric(10) ; x2[ c(5,6) ] = exp( c(-1,-0.25) )

oldpar = par( mfrow=c(2,2)  , omi=c(0,0,0,0), mai=c(0.45,0.5,0.1,0) )
mybarplot( x1 )   ; mybarplot( x2 )     #  row #1
mybarplot( 1-x1 ) ; mybarplot( 1-x2 )   #  row #2
par( oldpar )

## ----echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6.5, fig.height=4, fig.cap='Figure 2.2', out.width="100%", cache=FALSE----------------

mystepplot <- function( x )  {
# assumption: x is Type I
n = length(x)
plot( c(1/2,n+1/2), c(0,1), type='n', tcl=0, las=1, xlab='', ylab='', lab=c(n,5,7), mgp=c(3,0.25,0) )
grid( lty=1 )
beta = seq(1/2,n+1/2,by=1) ; segments( beta, 0, beta, -0.02 )
ij = which( 0<x & x<1)  ;  lambda = ij + c(1/2 - x[ ij[1] ], x[ ij[2] ] - 1/2)
lines( c(0.5,lambda[1]), c(0,0) ) ; lines(lambda,c(1,1)) ; lines( c(lambda[2],n+1/2), c(0,0) )
segments( lambda, c(0,0), lambda, c(1,1), lty=3 )
}

oldpar = par( mfrow=c(2,2), omi=c(0,0,0,0), mai=c(0.45,0.5,0.1,0) )
mystepplot( x1 ) ; mystepplot( x2 )     #  row #1
mybarplot( x1 ) ; mybarplot( x2 )       #  row #2
par( oldpar )

## ----rgl=TRUE, dev='png', echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6.5, fig.height=4, fig.cap='Figure 2.3  &emsp;&emsp;&emsp;  [these are interactive WebGL widgets]', fig.keep='none', fig.show='hide', out.width="100%", cache=FALSE----
rgl::par3d( zoom=0.7 )
rgl::mfrow3d( 1, 2 )
zono =  polarzonohedron(9)
plot2trans( zono )
rgl::next3d()
plot2trans( zono, level=c(0,4,7) )
rgl::rglwidget( webgl=TRUE )

## ----echo=FALSE,  message=TRUE,  warning=TRUE, fig.width=8, fig.height=3, fig.cap='Figure 10.1', out.width="100%", cache=FALSE----------------

plot_slabs <- function()
    {
    plot.new()

    xlim = c(-10,10)
    ylim = c(-7,7)

    theta   = 20 * pi/180
    rot2x2  = matrix( c(cos(theta),sin(theta),-sin(theta),cos(theta)), 2, 2 )

    plot.window( xlim, ylim, asp=1 )

    #   big slab
    x  = c(-15,15,15,-15)
    y   = c(-5,-5,5,5)
    xy  = rbind( x, y )

    xyrot   = rot2x2 %*% xy
    polygon( xyrot[1, ], xyrot[2, ], col='gray90' )

    xya = cbind( c(1,5), c(4,5) )
    xyrot   = rot2x2 %*% xya

    xymid   = rowMeans(xyrot)

    lines( xyrot[1, ], xyrot[2, ], lwd=5 )
    text( xymid[1], xymid[2], "abundant", adj=c(1,-1/2) )

    lines( -xyrot[1, ], -xyrot[2, ], lwd=5 )
    text( -xymid[1], -xymid[2], "abundant", adj=c(0,3/2) )

    #   small slab
    ytop    = 2.5
    x  = c(-15,15,15,-15)
    y   = c(-ytop,-ytop,ytop,ytop)
    xy  = rbind( x, y )

    xyrot   = rot2x2 %*% xy
    polygon( xyrot[1, ], xyrot[2, ], col='gray80', lty=2 )

    xyd = cbind( c(-ytop,ytop), c(0,ytop) )
    xyrot   = rot2x2 %*% xyd

    xymid   = rowMeans(xyrot)

    lines( xyrot[1, ], xyrot[2, ], lwd=5 )
    text( xymid[1], xymid[2], "deficient", adj=c(1,-1/2) )

    lines( -xyrot[1, ], -xyrot[2, ], lwd=5 )
    text( -xymid[1], -xymid[2], "deficient", adj=c(0,3/2) )

    xya     = cbind( c(-6,5), c(-6,5+2) )
    xyrot   = rot2x2 %*% xya

    arrows( xyrot[1,1], xyrot[2,1],  xyrot[1,2], xyrot[2,2], length=0.1, angle=20 )
    arrows( -xyrot[1,1], -xyrot[2,1],  -xyrot[1,2], -xyrot[2,2], length=0.1, angle=20 )
    
    #   label both slabs
    x0  = 7
    
    xy  = cbind( c( x0, (5+ytop)/2 ),  c( x0, -(5+ytop)/2 ) )
    xyrot   = rot2x2 %*% xy
    text( xyrot[1, ], xyrot[2, ], expression( S ) )

    xy  = c( x0, 0 )
    xyrot   = rot2x2 %*% xy
    text( xyrot[1], xyrot[2], expression( S[2] ) )

    points( 0, 0, pch=20 )

    #return( TRUE )
    }


plot_slab <- function()
    {
    plot.new()

    xlim = c(-10,10)
    ylim = c(-8,8)

    theta   = 20 * pi/180
    rot2x2  = matrix( c(cos(theta),sin(theta),-sin(theta),cos(theta)), 2, 2 )

    plot.window( xlim, ylim, asp=1 )

    #   big slab
    x  = c(-15,15,15,-15)
    y   = c(-5,-5,5,5)
    xy  = rbind( x, y )

    xyrot   = rot2x2 %*% xy
    polygon( xyrot[1, ], xyrot[2, ], col='gray80' )

    xya = cbind( c(1,5), c(4,5) )
    xyrot   = rot2x2 %*% xya

    xymid   = rowMeans(xyrot)

    lines( xyrot[1, ], xyrot[2, ], lwd=5 )
    text( xymid[1], xymid[2], "coincident", adj=c(1,-1/2) )

    lines( -xyrot[1, ], -xyrot[2, ], lwd=5 )
    text( -xymid[1], -xymid[2], "coincident", adj=c(0,3/2) )

    #   arrows
    xya     = cbind( c(-6,5), c(-6,5+2) )
    xyrot   = rot2x2 %*% xya

    arrows( xyrot[1,1], xyrot[2,1],  xyrot[1,2], xyrot[2,2], length=0.1, angle=20 )
    arrows( -xyrot[1,1], -xyrot[2,1],  -xyrot[1,2], -xyrot[2,2], length=0.1, angle=20 )

    
    #   label slab
    xy  = c( 6, 0 )
    xyrot   = rot2x2 %*% xy
    text( xyrot[1], xyrot[2], expression( S[2] == S ) )

    points( 0, 0, pch=20 )
    
    }
    
oldpar = par( mfrow=c(1,2)  , omi=c(0,0,0,0), mai=c(0,0.1,0,0.1) )

plot_slabs() ; plot_slab()

par( oldpar )

## ----echo=TRUE,  message=TRUE,  warning=TRUE--------------------------------------------------------------------------------------------------
matgen = colorimetry.genlist[[2]]   # the CIE 1931 CMFs at 1nm step
matgen = 100 * matgen / sum( matgen[2, ] )   # it's traditional to scale so the center has Y=50
zono =  zonohedron( matgen )
getcenter(zono) ; dim( getmatrix( getsimplified( getmatroid(zono) ) ) )
transitionsdf( zono )

## ----echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6.5, fig.height=3, fig.cap='Figure 10.2', out.width="100%", cache=FALSE---------------
oldpar = par( omi=c(0,0,0,0), mai=c(0.45,0.5,0.1,0) )
gnd = getground( getsimplified( getmatroid(zono) ) )
pcube = boundarypgramdata( zono, c(570,608), cube=TRUE )$pcube
xlim = range( gnd[which(0<pcube)] ) + 20*c(-1,1)
plot( xlim, c(0,1), type='n', xlab='', ylab='', las=1, lab=c(5,10,7), cex.axis=0.8 )
grid( col='gray', lty=1 )
lines( gnd, pcube, type='s' )
par( oldpar )

## ----rgl=TRUE, dev='png', echo=TRUE,  message=TRUE,  warning=FALSE, fig.width=6.5, fig.height=4, fig.cap='Figure 10.3', fig.keep='last', fig.show='hold', out.width="100%", cache=FALSE----
library( orientlib )
user3x3 = orientlib::rotmatrix( orientlib::eulerzyx( -0.249417, 0.7116067, 2.324364 ) )@x
dim(user3x3) = c(3,3)
par3d( userMatrix=rotationMatrix(matrix=user3x3), zoom=0.35 )
plothighertrans( zono )

## ----echo=FALSE, results='asis'-----------------------------------------------
options( old_opt )
sessionInfo()

