#include <R.h>
 #include <math.h>
 void p_1_shCTRL_x6pq2niv ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]*log(10.0)) ;
y[5+i**l] = exp(p[1]*log(10.0)) ;
y[7+i**l] = exp(p[2]*log(10.0)) ;
y[9+i**l] = 1.0 ;
y[10+i**l] = exp(p[3]*log(10.0)) ;
y[11+i**l] = exp(p[4]*log(10.0)) ;
y[13+i**l] = exp(p[5]*log(10.0)) ;
y[14+i**l] = exp(p[6]*log(10.0)) ;
y[15+i**l] = exp(p[7]*log(10.0)) ;
y[16+i**l] = exp(p[8]*log(10.0)) ;
y[17+i**l] = exp(p[9]*log(10.0)) ;
y[18+i**l] = exp(p[10]*log(10.0)) ;
y[19+i**l] = exp(p[11]*log(10.0)) ;
y[20+i**l] = p[12] ;
y[21+i**l] = exp(p[13]*log(10.0)) ;
y[22+i**l] = 1.0 ; 
}
}