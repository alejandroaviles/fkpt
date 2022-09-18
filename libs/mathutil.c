
#include "stdinc.h"
#include "mathfns.h"
#include "constant.h"
#include "mathutil.h"

void Rotation3D(real vec[], real alpha, real beta, real gamma)
{
    real a11, a12, a13;
    real a21, a22, a23;
    real a31, a32, a33;
	real vecp[3];
	real alp, bet, gam;
    int i, ndim;

	ndim=3;

	alp = radian(alpha);
	bet = radian(beta);
	gam = radian(gamma);

	a11 = rcos(gam)*rcos(bet)*rcos(alp) - rsin(gam)*rsin(alp);
	a12 = rcos(gam)*rcos(bet)*rsin(alp) + rsin(gam)*rcos(alp);
	a13 = -rcos(gam)*rsin(bet);

	a21 = -rsin(gam)*rcos(bet)*rcos(alp) - rcos(gam)*rsin(alp);
	a22 = -rsin(gam)*rcos(bet)*rsin(alp) + rcos(gam)*rcos(alp);
	a23 = rsin(gam)*rsin(bet);

	a31 = rsin(bet)*rcos(alp);
	a32 = rsin(bet)*rsin(alp);
	a33 = rcos(bet);

	vecp[0] = a11*vec[0] + a12*vec[1] + a13*vec[2];
	vecp[1] = a21*vec[0] + a22*vec[1] + a23*vec[2];
	vecp[2] = a31*vec[0] + a32*vec[1] + a33*vec[2];

    for (i = 0; i < ndim; i++)
        vec[i] = vecp[i];
}

real radian(real degree)
{
	return(degree * PI_D / ((real) 180.0) );
}

real degree(real radian)
{
	return(radian * ((real) 180.0) / PI_D );
}

real angle(real XI, real YI, real XF, real YF)
{
	real DX,DY, ang;

    DX = XF - XI;
    DY = YF - YI;
    if (DY == ZERO) {
		ang = ZERO;
        if ( DX < ZERO ) ang = ang + PI_D;
        return(ang);
    }
    if (DX == ZERO) {
		ang = PI_D/((real) 2.0);
        if ( DY < ZERO ) ang = ang + PI_D;
        return(ang);
	}
    ang = ratan( rabs(DY/DX) );
    if (DX < ZERO && DY > ZERO) {
        ang=PI_D-ang;
        return(ang);
	}
    if (DX < ZERO && DY < ZERO) {
		ang=PI_D+ang;
        return(ang);
	}
	if (DX > ZERO && DY < ZERO) {
        ang= ((real) 2.0) * PI_D - ang;
        return(ang);
	}
	return(ang);
}


double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}


double dmin(double x,double y)
{
  if(x<y)
    return x;
  else
    return y;
}


int imax(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}


int imin(int x,int y)
{
  if(x<y)
    return x;
  else
    return y;
}


void int_piksrt(int n, int arr[])	// Ordena un arreglo de enteros
{									// de menor a mayor.
	int i, j;
	int a;

	for (j=1; j<n; j++) {
		a = arr[j];
		i=j-1;
		while (i>=0 && arr[i] > a) {
			arr[i+1]=arr[i];
			i--;
		}
		arr[i+1]=a;
	}
}

// Sorting using the Heapsort Method (pr_anvorpol.c) (Ascending order)

void HeapSort (real *a, int *seq, int n)
{
  real q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] < a[seq[j + 1]]) ++ j;
        if (q < a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting an array of integers using the Heapsort Method (pr_anvorpol.c) (Ascending order)

void HeapSortInt (int *a, int *seq, int n)
{
  int q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] < a[seq[j + 1]]) ++ j;
        if (q < a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting using the Heapsort Method (Descending order)

void HeapSortDescend (real *a, int *seq, int n)
{
  real q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] > a[seq[j + 1]]) ++ j;
        if (q > a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting integer array using the Heapsort Method (Descending order)

void HeapSortIntDescend (int *a, int *seq, int n)
{
  int q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] > a[seq[j + 1]]) ++ j;
        if (q > a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Fast Fourier Transform (pr_anspcor.c)
// Use: FftComplex (a, n);
//  Cmplx *work;
//  AllocMem (work, 2 * (nValCorr - 1), Cmplx);
//  CSet (work[n], corrSum[j][k * nValCorr + n] * damp, 0.);
//  work[n] = work[2 * (nValCorr - 1) - n];
void FftComplex (Cmplx *a, int size)
{
  Cmplx t, w, wo;
  real theta;
  int i, j, k, n;

  k = 0;
  for (i = 0; i < size; i ++) {
    if (i < k) {
      t = a[i];
      a[i] = a[k];
      a[k] = t;
    }
    n = size / 2;
    while (n >= 1 && k >= n) {
      k -= n;
      n /= 2;
    }
    k += n;
  }
  for (n = 1; n < size; n *= 2) {
    theta = M_PI / n;
    CSet (wo, cos (theta) - 1., sin (theta));
    CSet (w, 1., 0.);
    for (k = 0; k < n; k ++) {
      for (i = k; i < size; i += 2 * n) {
        j = i + n;
        CMul (t, w, a[j]);
        CSub (a[j], a[i], t);
        CAdd (a[i], a[i], t);
      }
      CMul (t, w, wo);
      CAdd (w, w, t);
    }
  }
}

real Integrate (real *f, int nf)
{
  real s;
  int i;

  s = 0.5 * (f[0] + f[nf - 1]);
  for (i = 1; i < nf - 1; i ++) s += f[i];
  return (s);
}


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort3arrays(unsigned long n, double arr[], double brr[], double crr[])
{
    unsigned long i,ir=n,j,k,l=1,*istack;
    int jstack=0;
    double a,b,c,temp;
    
    istack=nr_lvector(1,NSTACK);
    for (;;) {
        if (ir-l < M) {
            for (j=l+1;j<=ir;j++) {
                a=arr[j];
                b=brr[j];
                c=crr[j];
                for (i=j-1;i>=l;i--) {
                    if (arr[i] <= a) break;
                    arr[i+1]=arr[i];
                    brr[i+1]=brr[i];
                    crr[i+1]=crr[i];
                }
                arr[i+1]=a;
                brr[i+1]=b;
                crr[i+1]=c;
            }
            if (!jstack) {
                free_lvector(istack,1,NSTACK);
                return;
            }
            ir=istack[jstack];
            l=istack[jstack-1];
            jstack -= 2;
        } else {
            k=(l+ir) >> 1;
            SWAP(arr[k],arr[l+1])
            SWAP(brr[k],brr[l+1])
            SWAP(crr[k],crr[l+1])
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir])
                SWAP(brr[l],brr[ir])
                SWAP(crr[l],crr[ir])
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir])
                SWAP(brr[l+1],brr[ir])
                SWAP(crr[l+1],crr[ir])
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l],arr[l+1])
                SWAP(brr[l],brr[l+1])
                SWAP(crr[l],crr[l+1])
            }
            i=l+1;
            j=ir;
            a=arr[l+1];
            b=brr[l+1];
            c=crr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j])
                SWAP(brr[i],brr[j])
                SWAP(crr[i],crr[j])
            }
            arr[l+1]=arr[j];
            arr[j]=a;
            brr[l+1]=brr[j];
            brr[j]=b;
            crr[l+1]=crr[j];
            crr[j]=c;
            jstack += 2;
            if (jstack > NSTACK) nrerror("NSTACK too small in sort3arrays.");
            if (ir-i+1 >= j-l) {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            } else {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
}
#undef M
#undef NSTACK
#undef SWAP


// Fitting

void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
         double *b, double *siga, double *sigb, double *chi2, double *q)
{
    double gammq(double a, double x);
    int i;
    double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
    
    *b=0.0;
    if (mwt) {
        ss=0.0;
        for (i=1;i<=ndata;i++) {
            wt=1.0/SQR(sig[i]);
            ss += wt;
            sx += x[i]*wt;
            sy += y[i]*wt;
        }
    } else {
        for (i=1;i<=ndata;i++) {
            sx += x[i];
            sy += y[i];
        }
        ss=ndata;
    }
    sxoss=sx/ss;
    if (mwt) {
        for (i=1;i<=ndata;i++) {
            t=(x[i]-sxoss)/sig[i];
            st2 += t*t;
            *b += t*y[i]/sig[i];
        }
    } else {
        for (i=1;i<=ndata;i++) {
            t=x[i]-sxoss;
            st2 += t*t;
            *b += t*y[i];
        }
    }
    *b /= st2;
    *a=(sy-sx*(*b))/ss;
    *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    *sigb=sqrt(1.0/st2);
    *chi2=0.0;
    *q=1.0;
    if (mwt == 0) {
        for (i=1;i<=ndata;i++)
            *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
        sigdat=sqrt((*chi2)/(ndata-2));
        *siga *= sigdat;
        *sigb *= sigdat;
    } else {
        for (i=1;i<=ndata;i++)
            *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
        if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(*chi2));
    }
}


