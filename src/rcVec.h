#include "CharSEXP.h"
#include "lessAndEqual.h" 



template <typename T>
class rcVec {		// a row vec or a col vec from a column-major order matrix
public:
    T * x; 		// pointer to the first element
    int len;    // length of vector: ncol for row vec; nrow for col vec
    int eltShift;  // index shift between adjacent elements: nrow for row vec; 1 for col vec
    int vecShift;  // index shift between adjacent vectors: 1 for row vec; nrow for col vec
    int nVec;		// number of vectors: nrow for row vec; ncol for col vec
    
    inline bool operator< (const rcVec& rhs ) const {
        // elementwise comparison of two vectors from the end
        // assuming equalTo<T>(usually operator==) and lessThan<T> (usually operator<) defined for type T 
        // also assuming operator= available for type T (Rcomplex is a struct of two doubles; SEXP in CharSEXP should be fine)
        T   L, R;
        for(int i=len-1; i>=0; i--){
            //  if ( equalTo<T>(L= *(x+eltShift*i), (R= *(rhs.x+rhs.eltShift*i))) ) continue;
            L   = *(x+eltShift*i);
            R   = *(rhs.x+rhs.eltShift*i);
            if( L == R )    continue ;
            return L < R ;   // lessThan<T>(L , R);  Glenn
        }
        return false;
    }
    inline bool operator== (const rcVec& rhs) const {
        T   L, R;        
        for(int i=len-1; i>=0; i--) {
            //  if ( ! equalTo<T>( *(x+eltShift*i),  *(rhs.x+rhs.eltShift*i)) ) return false;
            L   = *(x+eltShift*i);
            R   = *(rhs.x+rhs.eltShift*i);
            if( !(L == R) )    return false;    //  Glenn
        }
        return true;
	}
    /*
    friend inline bool operator> (const rcVec& lhs, const rcVec& rhs){return rhs < lhs;}
    friend inline bool operator<=(const rcVec& lhs, const rcVec& rhs){return !(lhs > rhs);}
    friend inline bool operator>=(const rcVec& lhs, const rcVec& rhs){return !(lhs < rhs);}
    */
};
