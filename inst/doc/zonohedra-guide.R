## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_opt = options( width=120 )

# if( !file.exists("figs") ) dir.create("figs")

require("rgl",quietly=TRUE)
rgl::setupKnitr(autoprint = TRUE)

## ----echo=TRUE,  message=FALSE----------------------------------------------------------------------------------------
library(zonohedra)
library(rgl)

## ----echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=8, fig.height=4, fig.cap='polar zonohedra with 5 generators (left) and 25 generators (right) &emsp;&emsp; [both of these are interactive WebGL widgets]', fig.keep='none', fig.show='hide', out.width="100%", cache=FALSE----
rgl::mfrow3d( 1, 2 )
pz5 = polarzonohedron( 5 ) ;  plot( pz5, ewd=5 )
rgl::next3d()
plot( polarzonohedron( 25 ), ewd=3 )
rgl::rglwidget( webgl=TRUE )

## ----echo=TRUE, message=FALSE-----------------------------------------------------------------------------------------
getmatrix( pz5 )

## ----echo=TRUE, message=FALSE-----------------------------------------------------------------------------------------
classics.genlist

## ----echo=TRUE, message=TRUE------------------------------------------------------------------------------------------
mat = classics.genlist[['TC']] ; mat

## ----rgl=TRUE, echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=8, fig.height=5, out.width="100%", fig.align="center", fig.cap='truncated cuboctahedron &emsp;&emsp;&emsp;&emsp; [This is an interactive WebGL widget]', fig.keep='last', fig.show='hide', cache=FALSE----
rgl::par3d( userMatrix = rotationMatrix( -20*pi/180, 0, 1, 1) )
zono = zonohedron( mat )
plot( zono, type='f' )
rgl::rglwidget( webgl=TRUE )

## ----echo=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------

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
pathtemp = tempdir()   #"./figs" ;   if( ! file.exists(pathtemp) ) dir.create(pathtemp)  # make temp folder
#  make a lot of .PNG files in pathtemp
movie3d( spin3d( getcenter(zono), rpm=revolutions*60/duration ), duration=duration, fps=fps, startTime=1/fps,
           convert=F, movie='junk', dir=pathtemp, verbose=F, webshot=F )
#  combine all the .PNGs into a single .webm
pathvec = dir( pathtemp, pattern="png$", full=T )

webm_file = sprintf( "%s/animation%g.webm", pathtemp, index )
out = av::av_encode_video( pathvec, output=webm_file, framerate=fps, codec='libvpx-vp9', verbose=F )
res = file.remove( pathvec )  # cleanup the .PNG files, leaving just the .webm

return( out )
}

## ----echo=TRUE, message=FALSE-----------------------------------------------------------------------------------------

video2html  <- function( path, attributes="controls loop autoplay muted" )
    {
    requireNamespace( "base64enc", quietly=TRUE )
    
    i   = regexpr( "[.][a-z]+$", path, ignore.case=T )
    if( i < 0 ) return('')
    ext = substring( path, i+1 )    # extract the extension, and skip over the '.'
    
    part1   = sprintf( '<video %s src="data:video/%s;base64,\n', attributes, ext )
    part2   = base64enc::base64encode( path, linewidth=120, newline='\n' )
    part3   = '"></video>'
    
    return( paste0( part1, part2, part3, collapse='' ) )
    }

## ----echo=TRUE, message=TRUE, warning=TRUE, fig.cap='object color solid', fig.keep='last', fig.show='hide', cache=FALSE----
# colorimetry.genlist[[1]] is a 3x81 matrix with the CIE 1931 CMFs at 5nm interval
zono5 = zonohedron( colorimetry.genlist[[1]] )
plot( zono5, type='f' )
webm_file = spinit( zono5, 2, vpsize=c(480,480) )

## ----echo=FALSE, message=TRUE, warning=TRUE---------------------------------------------------------------------------
video_html  = video2html(webm_file)
knitr::raw_html( video_html, meta=NULL, cacheable=FALSE )
unlink( dirname(webm_file) )

## ----echo=FALSE, results='asis'-----------------------------------------------
options( old_opt )
sessionInfo()

