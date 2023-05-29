#include <R.h>
#include <Rinternals.h>

//  #include <Defn.h>               // for LDOUBLE

#include <float.h>
#include <stdbool.h>

//  #include "config.h"         // this file is created by configure.win or configure

#include "macros.h"


/* Required by C99, but might be slow */
#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif


//  compute 3x3 determinant
double
det3x3( const double a[3], const double b[3], const double c[3] )
    {
    double  out = 0 ;

    out +=  a[0] * (b[1]*c[2] - b[2]*c[1]) ;
    out -=  b[0] * (a[1]*c[2] - a[2]*c[1]) ;
    out +=  c[0] * (a[1]*b[2] - a[2]*b[1]) ;

    return( out );
    }

//  a, b, c      points on the unit 2-sphere
//  for best numerics, the shortest side should be opposite c
//
//  returns *signed* area of the triangle.
//  positive area is counterclockwise, when viewed from the outside of the sphere

double
area_spherical_triangle( const double a[3], const double b[3], const double c[3] )
    {
    double  det = det3x3( a, b, c );

    //  Rprintf( "det = %g\n", det );

    if( fabs(det) < 5.e-16 )    return( 0.0 );  //  det too small

    double  sinhalfa    = 0;
    double  sinhalfb    = 0;
    double  sinhalfc    = 0;

    for( int k=0 ; k<3 ; k++ )
        {
        double  diffa, diffb, diffc;

        diffa   = b[k] - c[k] ;
        diffb   = a[k] - c[k] ;
        diffc   = a[k] - b[k] ;

        sinhalfa    += diffa*diffa ;
        sinhalfb    += diffb*diffb ;
        sinhalfc    += diffc*diffc ;
        }

    //  the next test is not needed; it is taken care of by the determinant test above
    //  if( sinhalfa==0 || sinhalfb==0 || sinhalfc==0 ) return(0.0);

    sinhalfa    = 0.5 * sqrt( sinhalfa );
    sinhalfb    = 0.5 * sqrt( sinhalfb );
    sinhalfc    = 0.5 * sqrt( sinhalfc );

    sinhalfa    = MIN2( sinhalfa, 1 );      //  necessity unknown, no known examples
    sinhalfb    = MIN2( sinhalfb, 1 );      //  necessity unknown, no known examples

    double  coshalfa    = sqrt( 1.0 - sinhalfa*sinhalfa );
    double  coshalfb    = sqrt( 1.0 - sinhalfb*sinhalfb );

    double  sinprod = sinhalfa * sinhalfb ;
    double  cosprod = coshalfa * coshalfb ;

    double  T2  = (sinhalfa*sinhalfa + sinhalfb*sinhalfb - 2*sinprod*sinprod - sinhalfc*sinhalfc) / (2*cosprod);

    double  cosC    = T2 / sinprod ;
    cosC            = MIN2( cosC, 1 );          //  necessary, testing uncovered an example
    double  sinC    = sqrt( 1.0 - cosC*cosC );

    double  E ;     //  spherical Excess

#if 0
    //  in the next line, the x part can be +, -, or 0.
    //  but the y part is always positive, and so E is between 0 and 2*M_PI, as it should be !
    E   = 2.0 * atan2( sinprod * sinC, cosprod + T2 ) ;        //  atan2f() is the same speed !
#else
    double  denom   = cosprod + T2;

    double  Ep ;
    if( denom != 0 )
        {
        Ep = 2.0 * atan( (sinprod * sinC) / denom ) ;

        if( Ep < 0 ) Ep += 2*M_PI;     // get a positive value using periodicity

        //if( 1.e-6 < fabs(E-Ep) )
        //    Rprintf( "y=%g  x=%g  E=%g   Ep=%g   E-Ep=%g\n", sinprod * sinC, cosprod + T2, E, Ep, E-Ep );
        }
    else
        {
        Ep = M_PI ; //  very special case
        }

    E   = Ep ;
#endif


#if 0
    if( E <= 0 )
        {
        Rprintf( "a=%g,%g,%g   b=%g,%g,%g   c=%g,%g,%g\n", a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
        Rprintf( "sinhalfa=%g  sinhalfb=%g  sinhalfc=%g\n", sinhalfa, sinhalfb, sinhalfc ) ;
        Rprintf( "E=%g < 0. sinprod=%g  cosC=%g   sinC=%g  cosprod=%g   T2=%g   temp=%g\n",
                        E, sinprod, cosC, sinC, cosprod, T2, temp ) ;
        }
#endif

    //  we want the *signed* area
    //  E   *= (0 < det) ? 1.0 : -1.0 ;
    if( det < 0 )   E = -E ;

    //  Rprintf( "signed area = %g\n", E );

    return( E ) ;
    }



SEXP
area_sphtri( SEXP sa, SEXP sb, SEXP sc )
    {
    SEXP    out = PROTECT( allocVector(REALSXP,1) );

    *(REAL(out)) = area_spherical_triangle( REAL(sa), REAL(sb), REAL(sc) ) ;

    UNPROTECT(1);

    return(out);
    }


//  load and unitize 4 rows of quadmat
bool
load_quadmat( const double offset[3], const double edge1[3], const double edge2[3],
                        double quadmat[4][3] )
    {
    int     j;

    for( j=0 ; j<3 ; j++ )
        {
        quadmat[0][j]   = offset[j] - 0.5 * edge1[j] - 0.5*edge2[j] ;
        quadmat[1][j]   = offset[j] - 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[2][j]   = offset[j] + 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[3][j]   = offset[j] + 0.5 * edge1[j] - 0.5*edge2[j] ;
        }

    //  normalize the columns of quadmat[][]
    //  project onto the unit sphere
    for( int i=0 ; i<4 ; i++ )
        {
        double  norm = 0 ;
        for( j=0 ; j<3 ; j++ )  norm += quadmat[i][j]*quadmat[i][j] ;

        if( fabs(norm) < 5.e-16 )
            return( false );

        norm    = sqrt(norm) ;

        for( j=0 ; j<3 ; j++ ) quadmat[i][j] /= norm;
        }

    return( true );
    }


//  smatgen     3xN matrix of generators, for the simplified matroid
//  sidxpair    N(N-1)/2 x 3 integer matrix of pairs of generators, 1-based
//  scenter     N(N-1)/2 x 3 matrix of pgram centers, in centered zonohedron coordinates
//  spoint      3D-point w.r.t. computing the linking number, must not be on a pgram, not verified
//                  this is in centered zonohedron coordinates
//
//  sidxpair and scenter do not have to be extended with antipodal data
//  this function does the extension automatically, whenever spoint is *not* 0 (the center of symmetry)

//  returns integer vector of length 1
//  method used: Gauss's integral definition

SEXP
linkingnumber( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP spoint )
    {
    const   int *dim ;

    dim = INTEGER(getAttrib(smatgen, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] < 3 )
        {
        Rprintf( "bad smatgen %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    int n = dim[1] ;
    const   double  *matgen = REAL(smatgen);

    int facets = (n*(n-1))/2 ;

    dim = INTEGER(getAttrib(sidxpair, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 2  )
        {
        Rprintf( "bad sidxpair %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   int *idxpair = INTEGER(sidxpair);

    dim = INTEGER(getAttrib(scenter, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 3  )
        {
        Rprintf( "bad scenter %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   double  *center = REAL(scenter);

    if( Rf_length(spoint) != 3 )    return(R_NilValue);
    const   double  *point = REAL(spoint);

    SEXP    out = PROTECT( allocVector(INTSXP,1) );
    *( INTEGER(out) ) = NA_INTEGER ;


    bool    symmetric = point[0]==0  &&   point[1]==0  &&  point[2]==0 ;


    //  quadmat[][] holds the vertices of the facet
    double  quadmat[4][3];

    LDOUBLE area = 0 ;

    //  Rprintf( "sizeof(area) = %d\n", sizeof(area) );

    for( int k=0 ; k<facets ; k++ )
        {
        double  offset[3];
        offset[0]   = center[k]             - point[0] ;
        offset[1]   = center[k + facets]    - point[1];
        offset[2]   = center[k + 2*facets]  - point[2] ;

        //  i and j are 1-based
        int i = idxpair[k] ;
        int j = idxpair[k + facets] ;

        const double    *edge1  = matgen + 3*(i-1) ;    //  1-based to 0-based
        const double    *edge2  = matgen + 3*(j-1) ;    //  1-based to 0-based

        if( ! load_quadmat(offset,edge1,edge2,quadmat) )
            {
            Rprintf( "linkingnumber(). The point (%g,%g,%g) (centered) is equal to a vertex of facet %d.\n",
                                point[0], point[1], point[2], k );
            Rprintf( "    The linking number is undefined; returning NA.\n" );
            UNPROTECT(1);
            return( out );
            }

        double  area_quads ;

        area_quads  = area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;

        //  Rprintf( "k=%d  area_quad=%g\n", k, area_quads );

#if 0
        if( k == 0 )
            {
            for( j=0 ; j<3 ; j++ )
                {
                Rprintf( "%10g %10g %10g %10g\n", quadmat[0][j],  quadmat[1][j], quadmat[2][j], quadmat[3][j] ) ;
                }
            Rprintf( "area[0] = %e\n", area );
            }
#endif

        if( ! symmetric )
            {
            //  the antipodal side will NOT give the same result
            //  repeat the calculations on the antipodal side
            //  on this pass, + point instead of -
            //  edge1 and edge2 stay the same
            offset[0]   = center[k]             + point[0] ;
            offset[1]   = center[k + facets]    + point[1] ;
            offset[2]   = center[k + 2*facets]  + point[2] ;

            if( ! load_quadmat(offset,edge1,edge2,quadmat) )
                {
                Rprintf( "linkingnumber(). The point (%g,%g,%g) (centered) is equal to a vertex of pgram %d.\n",
                                    point[0], point[1], point[2], k );
                Rprintf( "    The linking number is undefined; returning NA.\n" );
                UNPROTECT(1);
                return( out );
                }

            double  area_quad = area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;

            area_quads  += area_quad ;

            //  Rprintf( "k=%d  area_quad=%g\n", k, area_quad );
            }

        area += area_quads ;
        }

    if( symmetric )
        {
        //  the antipodal side will give the same result, so just double it
        area    *= 2.0 ;
        }

#if 0
    if( symmetric )
        {
        //  the antipodal side will give the same result, so just double it
        area *= 2.0 ;
        }
    else
        {
        //  the antipodal side will NOT give the same result
        //  repeat the calculations on the antipodal side
        for( int k=0 ; k<facets ; k++ )
            {
            double  offset[3];

            //  on this pass, + point instead of -
            offset[0]   = center[k]             + point[0] ;
            offset[1]   = center[k + facets]    + point[1] ;
            offset[2]   = center[k + 2*facets]  + point[2] ;

            //  i and j are 1-based
            int i = idxpair[k] ;
            int j = idxpair[k + facets] ;

            const double    *edge1  = matgen + 3*(i-1) ;    //  1-based to 0-based
            const double    *edge2  = matgen + 3*(j-1) ;    //  1-based to 0-based

            if( ! load_quadmat(offset,edge1,edge2,quadmat) )
                {
                Rprintf( "linkingnumber(). The point (%g,%g,%g) is equal to a vertex of facet %d.\n",
                                    point[0], point[1], point[2], k );
                UNPROTECT(1);
                return( out );
                }

            area += area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;
            }
        }
#endif

    //  4*M_PI is the area of the sphere
    //  change sign to be compatible with older conventions
    double  area_normalized = -area / (4*M_PI) ;

    int     linknum =  (int) roundf( area_normalized );

#if 1
    double  tol = 5.e-6 ;
    if( tol < fabs(area_normalized - linknum) )
        {
        Rprintf( "linkingnumber(). WARN.  fabs(area_normalized - linknum) = |%g|  >  %g (the tolerance).  Returning NA.\n",
                    area_normalized - linknum, tol );
        linknum = NA_INTEGER ;
        }
#endif

    *( INTEGER(out) ) = linknum ;

    UNPROTECT(1);

    return(out);
    }



//  smatcum     3 x N+1 matrix of the cumulative sum of the N generators, for the simplified matroid
//                  the first column is all 0s
//  spoint      3D-point w.r.t. computing the linking number, must not be on a pgram. not verified
//                  spoint is in *centered* zonohedron coordinates
//
//  returns integer vector of length 1
//  method used: Gauss's integral definition


SEXP
linkingnumber2( SEXP smatcum, SEXP spoint )
    {
    const   int *dim ;

    dim = INTEGER(getAttrib(smatcum, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] <= 3 )
        {
        Rprintf( "bad smatcum %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    int n = dim[1] - 1;
    const   double  *matcum = REAL(smatcum);

    //  verify that the first column is all 0s
    bool    ok = matcum[0]==0  &&  matcum[1]==0  &&  matcum[2]==0 ;
    if( ! ok )
        {
        Rprintf( "matcum is invalid; 1st column must be 0.\n" );
        return(R_NilValue);
        }

    if( Rf_length(spoint) != 3 )    return(R_NilValue);
    const   double  *point = REAL(spoint);

    const double    *white = matcum + 3*n ;

    double  center[3] ;
    for( int k=0 ; k<3 ; k++ )  center[k] = white[k]/2 ;

    //  allocate (N+1) x 2N x 3     3D array
    double  *vertex = Calloc( (n+1) * (2*n) * 3, double );

    int     rowstep = 3 ;
    int     colstep = 3*(n+1) ;

    //  fill vertex with the *centered* vertices.
    //  There are N(N-1) + 2 vertices in the surface, and N(N-1) + 2N cells in the array.

    //  start with "black" and "white", which are duplicated to many cells in the array
    //  this section fills 2N cells
    for( int jj=0 ; jj<n ; jj++ )
        {
        int     jb = 2*jj ;
        int     jw = jb + (n % 2) ;     // when n is odd, there is a shift

        double  *v0 = vertex  +  0*rowstep  +  jb*colstep ;
        double  *v1 = vertex  +  n*rowstep  +  jw*colstep ;

        for( int k=0 ; k<3 ; k++ )
            {
            v0[k]   = -center[k];            //  *black* is -center
            v1[k]   =  center[k];            //  *white* is +center
            }
        }

    int row, col;

    //  this section fills N(N-1) cells with *centered* vertices
    for( int j1=0 ; j1<n-1 ; j1++ )
        {
        const double    *s1 = matcum + j1*3;

        for( int j2=j1+1 ; j2<n ; j2++ )
            {
            const double    *s2 = matcum + j2*3;

            row = j2 - j1 ;
            col = j1 + j2 ;
            double  *v = vertex  +  row*rowstep  +  col*colstep ;

            //  now the antipodal cell
            row = n - row ;
            col = (col + n) % (2*n) ;
            double  *va = vertex  +  row*rowstep  +  col*colstep ;

            for( int k=0 ; k<3 ; k++ )
                {
                v[k]    = s2[k] - s1[k] - center[k] ;
                va[k]   = -v[k] ;
                }
            }
        }


    SEXP    out = PROTECT( allocVector(INTSXP,1) );
    *( INTEGER(out) ) = NA_INTEGER ;

    bool    symmetric = point[0]==0  &&  point[1]==0  &&  point[2]==0 ;

    //  in the next pass, subtract off point[] and unitize
    //  there are N(N-1) + 2N  cells

    int     rowmax = symmetric  ?  n/2+1 : n ;      // this opimization does not decrease time much - Jul 26 2022

    for( row=0 ; row<=rowmax ; row++ )
        {
        for( col=(row%2) ; col<2*n ; col+=2 )
            {
            double  *v = vertex  +  row*rowstep  +  col*colstep ;

            double  norm=0 ;
            for( int k=0 ; k<3 ; k++ )
                {
                v[k]    -= point[k] ;
                norm    += v[k]*v[k] ;
                }

            if( fabs(norm) < 5.e-16 )
                {
                Rprintf( "linkingnumber(). The point (%g,%g,%g) is equal to a vertex of the surface.\n",
                                    point[0], point[1], point[2] );
                Free( vertex );
                UNPROTECT(1);
                return( out );
                }

            norm    = sqrt(norm);

            for( int k=0 ; k<3 ; k++ )  v[k] /= norm ;
            }
        }


    //  iterate over the pgrams and compute area
    //  there are N(N-1) pgrams, unless symmetric when there are N(N-1)/2

    rowmax = symmetric  ?  n/2 : n-1 ;  // this opimization is effective - Jul 25 2022

    int     pgrams = 0 ;
    double  area = 0 ;
    for( row=1 ; row<=rowmax ; row++ )
        {
        int collim ;

        if( row<rowmax  ||  ! symmetric  ||  n%2 )
            collim  = 2*n ;
        else
            collim  = n ;

        for( col=((row+1)%2) ; col<collim ; col+=2 )
            {
            int colneg  = (col-1+2*n) % (2*n);
            int colpos  = (col+1) % (2*n);

            const double  *q0 = vertex  +  (row-1)*rowstep    +  col * colstep ;
            const double  *q1 = vertex  +  row*rowstep        +  colpos * colstep ;
            const double  *q2 = vertex  +  (row+1)*rowstep    +  col * colstep ;
            const double  *q3 = vertex  +  row*rowstep        +  colneg * colstep ;

            area += area_spherical_triangle( q1, q3, q0 )  +  area_spherical_triangle( q3, q1, q2 ) ;

            pgrams++ ;
            }
        }

    Free( vertex );

    if( symmetric )
        //  the antipodal side will give the same area sum, so just double it
        area    *= 2.0 ;


    //  4*M_PI is the area of the sphere
    //  change sign to be compatible with older conventions
    double  area_normalized = -area / (4*M_PI) ;

    int     linknum =  (int) roundf( area_normalized );

    *( INTEGER(out) ) = linknum ;

    UNPROTECT(1);       //  out


    //  Rprintf( "area_normalized=%g  area_normalized-linknum = %g\n", area_normalized, area_normalized-linknum );

#if 1
    int pgrams_corr = symmetric ? (n*(n-1))/2 : n*(n-1) ;

    if( pgrams != pgrams_corr )
        Rprintf( "ERROR. pgrams = %d  !=  %d (the correct value).\n", pgrams, pgrams_corr );

    double  tol = 5.e-7 ;
    if( tol < fabs(area_normalized - linknum) )
        Rprintf( "WARN. area_normalized - linknum = %g  >  %g\n", area_normalized - linknum, tol );
#endif

    return( out );
    }
