#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
(dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
(minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
(lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
(lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
(imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void sinft(float y[], int n);
void dsinft(double y[], int n);
void realft(float data[], unsigned long n, int isign);
void drealft(double data[], unsigned long n, int isign);
void four1(float data[], unsigned long nn, int isign);
void dfour1(double data[], unsigned long nn, int isign);

void gaussj(double **a, int n, double **b, int m);

void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb);
float zbrent(float (*func)(float), float x1, float x2, float tol);

float merff(float x);							
												
												
void nrerror(char error_text[]);
/* float *vector(long nl, long nh); */

float ran0(long *idum);
float fran1(long *idum);
double ran1(long *idum);
float ran2(long *idum);
float ran3(long *idum);
float ran4(long *idum);
float fgasdev(long *idum);
double gasdev(long *idum);

void psdes(unsigned long *lword, unsigned long *irword);

float *nr_vector(long nl, long nh);
int *nr_ivector(long nl, long nh);
unsigned char *nr_cvector(long nl, long nh);
unsigned long *nr_lvector(long nl, long nh);
/* real *nr_dvector(long nl, long nh); */
double *nr_dvector(long nl, long nh); 
double *dvector(long nl, long nh); 
float **nr_matrix(long nrl, long nrh, long ncl, long nch);
double **nr_dmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **nr_imatrix(long nrl, long nrh, long ncl, long nch);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

#endif /* _NR_UTILS_H_ */
