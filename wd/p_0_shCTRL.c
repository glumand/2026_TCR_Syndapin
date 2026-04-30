#include <R.h>
 #include <math.h>
 void p_0_shCTRL_71wxmsdc ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[7+i**l] = exp(p[0]*log(10.0)) ;
y[10+i**l] = exp(p[1]*log(10.0)) ;
y[11+i**l] = exp(p[2]*log(10.0)) ;
y[12+i**l] = exp(p[3]*log(10.0)) ;
y[13+i**l] = exp(p[4]*log(10.0)) ;
y[14+i**l] = exp(p[5]*log(10.0)) ;
y[15+i**l] = exp(p[6]*log(10.0)) ;
y[16+i**l] = 1.0 ;
y[17+i**l] = exp(p[7]*log(10.0)) ;
y[18+i**l] = p[8] ; 
}
}