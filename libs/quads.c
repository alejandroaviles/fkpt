

#include "stdinc.h"
#include "mathfns.h"
//#include "../general/constant.h"
#include "quads.h"
#include "numrec.h"

#define EPS 3.0e-11

void gauleg(double x1, double x2, double x[], double w[], int n)
{
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        z=cos(3.141592654*(i-0.25)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }
}
#undef EPS

double qgauss(double (*func)(double), double a, double b,
              double x[], double w[], int n)
{
    int j;
    double s;
    
    s=0;
    for (j=1;j<=n;j++) {
        s += w[j]*(*func)(x[j]);
    }
    return s;
}


// qromo

//#include <math.h>
//#define EPS 1.0e-6 // Original
#define JMAX 14
#define JMAXP (JMAX+1)
//#define K 5

//double qromo(double (*func)(double), double a, double b,
//             double (*choose)(double(*)(double), double, double, int), double epsq)
//double (*choose)(double(*)(double), double, double, int)) // Original
double qromo(double (*func)(double), double a, double b,
             double (*choose)(double(*)(double), double, double, int),
             double epsq, int KK)
{
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    void nrerror(char error_text[]);
    int j;
    double ss,dss,h[JMAXP+1],s[JMAXP];
    
    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
        s[j]=(*choose)(func,a,b,j);
        if (j >= KK) {
            polint(&h[j-KK],&s[j-KK],KK,0.0,&ss,&dss);
            //            if (fabs(dss) <= EPS*fabs(ss)) return ss; // Original
            if (fabs(dss) <= epsq*fabs(ss)) return ss;
        }
        h[j+1]=h[j]/9.0;
    }
    nrerror("Too many steps in routing qromo");
    return 0.0;
}
//#undef EPS // Original
#undef JMAX
#undef JMAXP
//#undef K

// END qromo


// polint


void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;
    double *c,*d;
    
    dif=fabs(x-xa[1]);
    c=dvector(1,n);
    d=dvector(1,n);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free_dvector(d,1,n);
    free_dvector(c,1,n);
}

// END polint



// BEGIN midexp

#define FUNC(x) ((*funk)(-log(x))/(x))

double midexp(double (*funk)(double), double aa, double bb, int n)
{
    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;
    
    b=exp(-aa);
    a=0.0;
    if (n == 1) {
        return (s=(b-a)*FUNC(0.5*(a+b)));
    } else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}
#undef FUNC
// END midexp



// BEGIN midinf

#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))

double midinf(double (*funk)(double), double aa, double bb, int n)
{
    double x,tnm,sum,del,ddel,b,a;
    static double s;
    int it,j;
    
    b=1.0/aa;
    a=1.0/bb;
    if (n == 1) {
        return (s=(b-a)*FUNC(0.5*(a+b)));
    } else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        return (s=(s+(b-a)*sum/tnm)/3.0);
    }
}
#undef FUNC
// END midinf



// BEGIN midsql

//#include <math.h>
#define FUNC(x) (2.0*(x)*(*funk)(aa+(x)*(x)))

double midsql(double (*funk)(double), double aa, double bb, int n)
{
    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;
    
    b=rsqrt(bb-aa);
    a=0.0;
    if (n == 1) {
        return (s=(b-a)*FUNC(0.5*(a+b)));
    } else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}
#undef FUNC
// END midsql



// BEGIN midsqu

//#include <math.h>
#define FUNC(x) (2.0*(x)*(*funk)(bb-(x)*(x)))

double midsqu(double (*funk)(double), double aa, double bb, int n)
{
    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;
    
    b=rsqrt(bb-aa);
    a=0.0;
    if (n == 1) {
        return (s=(b-a)*FUNC(0.5*(a+b)));
    } else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}
#undef FUNC
// END midsqu

// END midpnt
#define FUNC(x) ((*func)(x))

double midpnt(double (*func)(double), double a, double b, int n)
{
    double x,tnm,sum,del,ddel;
    static double s;
    int it,j;
    
    if (n == 1) {
        return (s=(b-a)*FUNC(0.5*(a+b)));
    } else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}
#undef FUNC
// END midpnt


void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
    int i,k;
    double p,qn,sig,un,*u;
    
    u=dvector(1,n-1);
    if (yp1 > 0.99e30)
        y2[1]=u[1]=0.0;
    else {
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    free_dvector(u,1,n-1);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
    void nrerror(char error_text[]);
    int klo,khi,k;
    double h,b,a;
    
    klo=1;
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void locate(double xx[], unsigned long n, double x, unsigned long *j)
{
    unsigned long ju,jm,jl;
    int ascnd;

    jl=0;
    ju=n+1;
    ascnd=(xx[n] >= xx[1]);
    while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    if (x == xx[1]) *j=1;
    else if(x == xx[n]) *j=n-1;
    else *j=jl;
}
