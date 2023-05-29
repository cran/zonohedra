## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_opt = options( width=120 )

require("rgl",quietly=TRUE)
rgl::setupKnitr(autoprint = TRUE)

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------
library(zonohedra)

## ---- echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6, fig.height=4, fig.cap='', out.width="80%", cache=FALSE----
zono =  polarzonogon( 14, 4 )
oldpar = par( omi=c(0,0,0,0), mai=c(0.8,0.7,0.7,0.2) )
plot( zono, elabels=T )
par( oldpar )

## ---- echo=TRUE,  message=TRUE,  warning=TRUE, fig.width=6, fig.height=4, fig.cap='', out.width="80%", cache=FALSE----
oldpar = par( omi=c(0,0,0,0), mai=c(0.8,0.7,0.7,0.2) )
plot( zono, tiling=T, elabels=T, tlabels=T )
par( oldpar )

## ---- echo=FALSE, results='asis'----------------------------------------------
options( old_opt )
sessionInfo()

