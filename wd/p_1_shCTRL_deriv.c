#include <R.h>
 #include <math.h>
 void p_1_shCTRL_deriv_f098ph95 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[1+i**l] = exp(p[0]*log(10.0))*log(10.0) ;
y[22+i**l] = exp(p[1]*log(10.0))*log(10.0) ;
y[47+i**l] = exp(p[2]*log(10.0))*log(10.0) ;
y[70+i**l] = exp(p[3]*log(10.0))*log(10.0) ;
y[91+i**l] = exp(p[4]*log(10.0))*log(10.0) ;
y[112+i**l] = exp(p[5]*log(10.0))*log(10.0) ;
y[133+i**l] = exp(p[6]*log(10.0))*log(10.0) ;
y[154+i**l] = exp(p[7]*log(10.0))*log(10.0) ;
y[175+i**l] = exp(p[8]*log(10.0))*log(10.0) ;
y[197+i**l] = exp(p[9]*log(10.0))*log(10.0) ;
y[218+i**l] = 1.0 ; 
}
}