## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_opt = options( width=120 )

# if( !file.exists("figs") ) dir.create("figs")

require("rgl",quietly=TRUE)
rgl::setupKnitr(autoprint = TRUE)

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------
library(zonohedra)
library(rgl)

## ---- echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=8, fig.height=4, fig.cap='polar zonohedra with 5 generators (left) and 25 generators (right) &emsp;&emsp; [both of these are interactive WebGL widgets]', fig.keep='none', fig.show='hide', out.width="100%", cache=FALSE----
rgl::mfrow3d( 1, 2 )
pz5 = polarzonohedron( 5 ) ;  plot( pz5, ewd=5 )
rgl::next3d()
plot( polarzonohedron( 25 ), ewd=3 )
rgl::rglwidget( webgl=TRUE )

## ---- echo=TRUE, message=FALSE----------------------------------------------------------------------------------------
getmatrix( pz5 )

## ---- echo=TRUE, message=FALSE----------------------------------------------------------------------------------------
classics.genlist

## ---- echo=TRUE, message=TRUE-----------------------------------------------------------------------------------------
mat = classics.genlist[['TC']] ; mat

## ---- rgl=TRUE, echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=8, fig.height=5, out.width="100%", fig.align="center", fig.cap='truncated cuboctahedron &emsp;&emsp;&emsp;&emsp; [This is an interactive WebGL widget]', fig.keep='last', fig.show='hide', cache=FALSE----
rgl::par3d( userMatrix = rotationMatrix( -20*pi/180, 0, 1, 1) )
zono = zonohedron( mat )
plot( zono, type='f' )
rgl::rglwidget( webgl=TRUE )

## ---- echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------
library(gifski)

#   zono        the zonohedron
#   id          unique ID for this animation, a positive integer
#   fps         frames per second
#   duration    of the animation, in seconds
#   revolutions number of revolutions
#   vpsize      viewport size = (width,height)
spinit <- function( zono, index, fps=5, duration=8, revolutions=1, vpsize=c(480,480) ) {
#  enlarge viewport
wr = par3d( "windowRect" ) 
par3d( windowRect = c( wr[1:2], wr[1:2] + vpsize ) )
pathtemp = "./figs" ;   if( ! file.exists(pathtemp) ) dir.create(pathtemp)  # make temp folder
#  make a lot of .PNG files in pathtemp
movie3d( spin3d( getcenter(zono), rpm=revolutions*60/duration ), duration=duration, fps=fps, startTime=1/fps,
           convert=F, movie='junk', dir=pathtemp, verbose=F, webshot=F )
#  combine all the .PNGs into a single .GIF
pathvec = dir( pathtemp, pattern="png$", full=T )
gif_file = sprintf( "./figs/animation%g.gif", index ) 
# if( file.exists(gif_file) )  file.remove( gif_file )
out = gifski::gifski( pathvec, gif_file=gif_file, delay=1/fps, progress=F, width=vpsize[1], height=vpsize[2] )
res = file.remove( pathvec )  # cleanup the .PNG files, leaving just the .GIF

return( out )
}

## ---- echo=TRUE, message=TRUE, warning=TRUE, fig.cap='optimal color solid', fig.keep='last', fig.show='hide', cache=FALSE----
# colorimetry.genlist[[1]] is a 3x81 matrix with the CIE 1931 CMFs at 5nm interval
zono5 = zonohedron( colorimetry.genlist[[1]] )
plot( zono5, type='f' )
gif_file = spinit( zono5, 2, vpsize=c(256,256) )

## ---- echo=FALSE, results='asis'----------------------------------------------
options( old_opt )
sessionInfo()

