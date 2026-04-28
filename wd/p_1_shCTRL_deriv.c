#include <R.h>
 #include <math.h>
 void p_1_shCTRL_deriv_xul5a48a ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]*log(10.0))*log(10.0) ;
y[31+i**l] = exp(p[1]*log(10.0))*log(10.0) ;
y[59+i**l] = exp(p[2]*log(10.0))*log(10.0) ;
y[88+i**l] = exp(p[3]*log(10.0))*log(10.0) ;
y[115+i**l] = exp(p[4]*log(10.0))*log(10.0) ;
y[143+i**l] = exp(p[5]*log(10.0))*log(10.0) ;
y[170+i**l] = exp(p[6]*log(10.0))*log(10.0) ;
y[197+i**l] = exp(p[7]*log(10.0))*log(10.0) ;
y[224+i**l] = exp(p[8]*log(10.0))*log(10.0) ;
y[251+i**l] = exp(p[9]*log(10.0))*log(10.0) ;
y[278+i**l] = exp(p[10]*log(10.0))*log(10.0) ;
y[305+i**l] = exp(p[11]*log(10.0))*log(10.0) ;
y[332+i**l] = 1.0 ;
y[359+i**l] = exp(p[13]*log(10.0))*log(10.0) ; 
}
}