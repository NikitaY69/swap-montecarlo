#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include "nrutil.h"
#define ITMAX 10000
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,0,n);free_dvector(h,0,n);free_dvector(g,0,n);
#define TOL 2.0e-4
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SignR(x,y) ( ((y)>0)  ? (x) : (-x))
#define Sqr(x) ((x) * (x))
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define N 2
#define CUT 1.25
#define C0  -1.924145348608
#define C1  2.111062325330
#define C2  -0.591097451092
#define NSOFT 12
#define EPSP  0.2
#define SIZE 20000 // 4000
#define SIDE 100.0 //44.721359550000003


int ncom;
double *pcom,*xicom,(*nrfunc)(double []);
void (*nrdfun)(double [], double []);
extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);
extern void (*nrdfun)(double [], double []);
double sigma[SIZE/N] ;
double Pshift(double a);

typedef struct{
  double sigmaB;
  int **list;
  int  *nlist;
}strData;

strData data;

double Pshift(double a){
    return a - SIDE*floor((a+SIDE/2)/SIDE);
}

double definefunc(double p[]){
  int i,i2, j, k, n=SIZE, n3=SIZE/N;
  double fret=0, dist, rij, sigmaij;
  double half_side=SIDE*0.5;
  
  for( i=0; i<n3; i++){
    //for( j=i+1; j<n3; j++){
    for(i2=0; i2<data.nlist[i]; i2++){
      j=data.list[i][i2];
      sigmaij = (sigma[i]+sigma[j])*0.5* (1-EPSP*fabs(sigma[i]-sigma[j]));
      rij=0;
      for(k=0; k<N; k++){
        // p[i*N+k] = Pshift(p[i*N+k]); 
        // p[j*N+k] = Pshift(p[j*N+k]);
	dist = Pshift(p[i*N+k]) - Pshift(p[j*N+k]);
	if( fabs(dist) > half_side )
	  dist -= SignR(SIDE,dist);
	rij+= dist*dist;
      }
      
      rij=sqrt(rij);
      rij/=sigmaij;
      if(rij<CUT){
	fret+= pow(rij,-NSOFT) + C0 + C1*pow(rij,2)  + C2*pow(rij,4)  ;
      }
    }
  }
  return fret/((double)n3*2);
}

void definedfunc(double p[], double xi[]){
  int j, i, i2, k, n=SIZE, n3=SIZE/N;
  double half_side=SIDE*0.5;
  double dist, rij, sigmaij, rs, rdiff[N];

  for(i=0; i<n3; i++){
    xi[N*i]=0;      xi[N*i+1]=0;      // xi[N*i+2]=0;
    //for(j=0; j<n3; j++){
    for(i2=0; i2<data.nlist[i]; i2++){
      j=data.list[i][i2];
      if(j!=i){
	sigmaij = (sigma[i]+sigma[j])*0.5 *(1-EPSP*fabs(sigma[i]-sigma[j]));
	rij=0;
	for(k=0; k<N; k++){
    // p[i*N+k] = Pshift(p[i*N+k]); 
    // p[j*N+k] = Pshift(p[j*N+k]);
	  dist = Pshift(p[i*N+k]) - Pshift(p[j*N+k]); // ri-rj (xi,yi,zi)
	  if( fabs(dist) > half_side )
	    dist -= SignR(SIDE,dist);
	  rij+= dist*dist;
	  rdiff[k]=dist;
	}
	rij=sqrt(rij);
	rs=rij/sigmaij;
	if(rs<CUT){
	  for(k=0; k<N; k++){
	    xi[N*i+k] += (-NSOFT * pow(rs,-NSOFT) * pow(rij,-2) + 2*C1*pow(sigmaij, -2) + 4*C2*pow(sigmaij, -4)*rij*rij)*rdiff[k];
	  } // analycal gradient
	}
      }
    }
    for(k=0; k<N; k++){
      xi[N*i+k] /= (2.*n3);
    }
  }
}



void createList(double p[]){
  int i, j, k, part;
  double dist, dist2, side, sigmaAB, cut, cutPart, half_side, skin;

  part = SIZE/N;
  side = SIDE;
  half_side = SIDE*0.5;
  cut = CUT;
  skin = CUT*0.5;


  for(i=0; i<part; i++){
    data.nlist[i]=0;
  }
  for(i=0; i<part; i++){  
    sigmaAB = data.sigmaB;
    cutPart = Sqr(cut*sigmaAB + skin);

    for(j=i+1; j<part; j++){
      dist2=0;
      for(k=0;  k<N;   k++){
	dist = Pshift(p[i*N+k]) - Pshift(p[j*N+k]);
        if( fabs(dist) > half_side )
          dist -= SignR(side,dist);
        dist2 += Sqr(dist);
      }
      if( dist2 < cutPart){
	      data.list[i][data.nlist[i]]=j;
        data.list[j][data.nlist[j]]=i;
        data.nlist[i]++;
        data.nlist[j]++;
      }
    }
  }
}




int main(int argc, char *argv[]){
  void frprmn(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []));
  void readConf(double p[], char *input);
  void writeConf(double p[], char *output);
  int n=SIZE, n3=SIZE/N, iter, i;
  double p[n], fret;
  double ftol, initial;
  double (*func)(double []);
  void (*dfunc)(double [], double []);

  //  if ( (data = (strData *)  calloc(1, sizeof(strData))  ) == NULL ) exit(1);
  data.sigmaB=1.62;
  if( (data.nlist=(int*) calloc((SIZE/N),sizeof(int )))==NULL)exit(1);
  if(   (data.list   =(int**) calloc((SIZE/N),sizeof(int*)))==NULL)exit(1);
  for(i=0; i<n3; i++){
    //if( (data.list[i]=(short int*)  calloc((SIZE/N),sizeof(short int)))==NULL)exit(1);
    if( (data.list[i]=(int*)  calloc((600),sizeof(int)))==NULL)exit(1);
  }
    
  
  func= &definefunc;
  dfunc= &definedfunc;
  char *input = argv[1];
  char *output = argv[2];
  printf("IS     EQ   iterations\n");
    readConf(p, input);
    createList(p);
    initial = func(p);
    ftol=0.0000000001;
    frprmn(p,n, ftol, &iter, &fret, func, dfunc);
    //    printf("%i   %.12e   %.12e  %i\n", i, fret/(double)(n3), initial/(double)(n3), iter );
    printf("%.12e   %.12e  %i\n", fret, initial, iter );
    writeConf(p, output); 
  }


void frprmn(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double [])){
  void linmin(double p[], double xi[], int n, double *fret,double (*func)(double []) );
  void dlinmin(double p[], double xi[], int n, double *fret,double (*func)(double []), void (*dfunc)(double [], double []));
  int j,its, i, n3=SIZE/N;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;
  
  g=dvector(0,n-1);
  h=dvector(0,n-1);
  xi=dvector(0,n-1);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  //printf("#iteration 0  EQene  %.12e    %.12e\n", fp, *fret);
  for (j=0;j<n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    //linmin(p,xi,n,fret,func);
    dlinmin(p,xi,n,fret,func, dfunc);
    printf("iteration %i  value  %e    %e \n", its,  fp, *fret);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      FREEALL
	return;
    }

    fp= *fret;
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=0;j<n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      FREEALL
	return;
    }
    gam=dgg/gg;
    for (j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
    /*   for(i=0; i<n3; i++){
    //      printf("iter  %i   conf %f  %f  %f      %f  %f  %f\n", its,  p[i*N], p[i*N+1], p[i*N+2], xi[i*N], xi[i*N+1], xi[i*N+2]);
      }*/
    //printf("iteration %i  value  %f    %f    %f\n", its,  fp, *fret, func(p));
    fflush(0);
  }
  nrerror("Too many iterations in frprmn");
}






void linmin(double p[], double xi[], int n, double *fret, double (*func)(double [])){
  double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
  double f1dim(double x);
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,double *fc, double (*func)(double));
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n;
  pcom=dvector(0,n-1);
  xicom=dvector(0,n-1);
  nrfunc=func;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }

  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);

  //  printf("linmin  fa  %e   fx  %e  fb  %e\n", fa, fx,  fb);

  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);

  //  printf("linmin  fret  %e \n", *fret);

  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
    // if(p[j]>0.5*SIDE){
    //   p[j]=p[j]-SIDE;
    // }
    // if(p[j]<-0.5*SIDE){
    //   p[j]=SIDE-p[j];
    // }
  }
  free_dvector(xicom,0,n-1);
  free_dvector(pcom,0,n-1);
}


double f1dim(double x){
  int j;
  double f,*xt;
  xt=dvector(0,ncom-1);
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  //  printf("f1dim  f  %e  x  %e   %e  %e  %e\n", f, x, xt[0], xt[1], xt[2]);
  free_dvector(xt,0,ncom-1);
  return f;
}


double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin){
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  nrerror("Too many iterations in brent");
  *xmin=x;
  return fx;
}



void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double)){
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  //  printf("mnbrack00  fa  %e   fb  %e  fc  %e   fu  %e\n", *fa, *fb,  *fc, fu);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  //  printf("mnbrack0  fa  %e   fb  %e  fc  %e   fu  %e\n", *fa, *fb,  *fc, fu);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    //    printf("mnbrack1  fa  %e   fb  %e  fc  %e   fu  %e\n", *fa, *fb,  *fc, fu);
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
    //    printf("mnbrack2  fa  %e   fb  %e  fc  %e   fu  %e\n", *fa, *fb,  *fc, fu);
  }
  //  printf("mnbrack3  fa  %e   fb  %e  fc  %e   fu  %e\n", *fa, *fb,  *fc, fu);
}



void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []), void (*dfunc)(double [], double [])){
  double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double), double tol, double *xmin);
  double f1dim(double x);
  double df1dim(double x);
  void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n;
  pcom=dvector(0,n-1);
  xicom=dvector(0,n-1);
  nrfunc=func;
  nrdfun=dfunc;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
    // if(p[j]>0.5*SIDE){
    //   p[j]=p[j]-SIDE;
    // }
    // if(p[j]<-0.5*SIDE){
    //   p[j]=SIDE+p[j];
    // }
  }
  free_dvector(xicom,0,n);
  free_dvector(pcom,0,n);
}


double df1dim(double x){
  int j;
  double df1=0.0;
  double *xt,*df;
  xt=dvector(0,ncom-1);
  df=dvector(0,ncom-1);
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  (*nrdfun)(xt,df);
  for (j=0;j<ncom;j++) df1 += df[j]*xicom[j];
  free_dvector(df,0,ncom-1);
  free_dvector(xt,0,ncom-1);
  return df1;
}





double dbrent(double ax, double bx, double cx, double (*f)(double),double (*df)(double), double tol, double *xmin){
  int iter,ok1,ok2;
  double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1=2.0*(b-a);
      d2=d1;
      if (dw != dx) d1=(w-x)*dx/(dx-dw);
      if (dv != dx) d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e;
      e=d;
      if (ok1 || ok2) {
	if (ok1 && ok2)
	  d=(fabs(d1) < fabs(d2) ? d1 : d2); 
	else if (ok1)
	  d=d1;
	else
	  d=d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=SIGN(tol1,xm-x);
	} else {
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SIGN(tol1,d);
      fu=(*f)(u);
      if (fu > fx) {
	*xmin=x;
	return fx;
      }
    }
    du=(*df)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw)
	MOV3(w,fw,dw, x,fx,dx)
	MOV3(x,fx,dx, u,fu,du)
	} else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw)
	  MOV3(w,fw,dw, u,fu,du)
	  } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du)
	  }
    }
  }
  nrerror("Too many iterations in routine dbrent");
  return 0.0;
}




 void readConf(double p[], char *input){
  int i, j, lego, n=SIZE, n3=SIZE/N;
  double vacuum;
  FILE *fp;

  //printf("%s\n", filename);
  if(  (fp = fopen(input, "r"))==NULL){ printf("Errore apertura file configurazione\n"); exit(1);}
  lego = 0;
  

  for(j=0; j<n3; j++){
    lego+=fscanf(fp,"%lf",  &sigma[j]);
    for(i=0; i<N; i++){
      lego+=fscanf(fp,"%lf", &p[N*j+i]);
      // p[N*j+i] = Pshift(p[N*j+i]);
    };
    //for(i=0; i<N; i++){
     // lego+=fscanf(fp,"%lf", &vacuum);
   // }
    //    p[n+j] = sigma[j];
  }
  if(lego != (N+1) * n3){ printf("Errore lettura file configurazioni  lego:%i\n", lego); exit(1);}
  fclose(fp);  
}



 void writeConf(double p[], char *output){
  int i, j, lego, n=SIZE, n3=SIZE/N;
  double vacuum;
  FILE *fp;

  //printf("%s\n", filename);
  if(  (fp = fopen(output, "w"))==NULL){ printf("Errore apertura file configurazione\n"); exit(1);}
  lego = 0;
  
  // lego += fprintf(fp,"%i\n", SIZE/N);
  // lego += fprintf(fp,"%.12lf\n", SIDE );
  
  for(j=0; j<n3; j++){
    lego+=fprintf(fp,"%.12lf ",  sigma[j]);
    for(i=0; i<N; i++){
      lego+=fprintf(fp,"%.12lf ", p[N*j+i]);
    } lego+=fprintf(fp,"\n");
  }
  fclose(fp);  
}
