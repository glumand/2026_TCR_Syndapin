#include <R.h>
 #include <math.h>
 void p_0_shCTRL_deriv_ckbhhnec ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[7+i**l] = exp(p[0]*log(10.0))*log(10.0) ;
y[30+i**l] = exp(p[1]*log(10.0))*log(10.0) ;
y[51+i**l] = exp(p[2]*log(10.0))*log(10.0) ;
y[72+i**l] = exp(p[3]*log(10.0))*log(10.0) ;
y[93+i**l] = exp(p[4]*log(10.0))*log(10.0) ;
y[114+i**l] = exp(p[5]*log(10.0))*log(10.0) ;
y[135+i**l] = exp(p[6]*log(10.0))*log(10.0) ;
y[157+i**l] = exp(p[7]*log(10.0))*log(10.0) ;
y[178+i**l] = 1.0 ; 
}
}