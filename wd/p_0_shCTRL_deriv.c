#include <R.h>
 #include <math.h>
 void p_0_shCTRL_deriv_ms4oawpa ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[7+i**l] = exp(p[0]*log(10.0))*log(10.0) ;
y[36+i**l] = exp(p[1]*log(10.0))*log(10.0) ;
y[63+i**l] = exp(p[2]*log(10.0))*log(10.0) ;
y[91+i**l] = exp(p[3]*log(10.0))*log(10.0) ;
y[118+i**l] = exp(p[4]*log(10.0))*log(10.0) ;
y[145+i**l] = exp(p[5]*log(10.0))*log(10.0) ;
y[172+i**l] = exp(p[6]*log(10.0))*log(10.0) ;
y[199+i**l] = exp(p[7]*log(10.0))*log(10.0) ;
y[226+i**l] = exp(p[8]*log(10.0))*log(10.0) ;
y[253+i**l] = exp(p[9]*log(10.0))*log(10.0) ;
y[280+i**l] = 1.0 ;
y[307+i**l] = exp(p[11]*log(10.0))*log(10.0) ; 
}
}