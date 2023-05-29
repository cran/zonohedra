## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_opt = options( width=144 )

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------------------------------
library(zonohedra)

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------------------------------
matgen = colorimetry.genlist[[2]]   # the CIE 1931 CMFs at 1nm step
matgen = 100 * matgen / sum( matgen[2, ] )   # it is traditional to scale so the center has Y=50, recall we use Illuminant E
zono =  zonohedron( matgen )
base = getcenter(zono) ; base

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------------------------------
theta = 1.478858 ; phi = 0.371322
u = c( sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi) ) ; u

## ---- echo=TRUE,  message=TRUE----------------------------------------------------------------------------------------------------------------
df_opt = raytrace( zono, base, u ) ; df_opt
xyz_opt = df_opt$point[1, ] ; xyz_opt

## ---- echo=TRUE,  message=TRUE----------------------------------------------------------------------------------------------------------------
invertboundary( zono, xyz_opt )$transitions

## ---- echo=TRUE,  message=TRUE----------------------------------------------------------------------------------------------------------------
df_2trans = raytrace2trans( zono, base, u ) ; df_2trans
xyz_2trans = df_2trans$point[1, ] ; xyz_2trans

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------------------------------
df_opt$tmax - df_2trans$tmax

## ---- echo=TRUE,  message=FALSE---------------------------------------------------------------------------------------------------------------
xyz_mid = (xyz_opt + xyz_2trans) / 2
inside( zono, xyz_mid )
inside2trans( zono, xyz_mid )

## ---- echo=FALSE, results='asis'----------------------------------------------
options(old_opt)
sessionInfo()

