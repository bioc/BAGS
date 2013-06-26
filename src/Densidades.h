/*
 *  Densidades.h
 *
 *  Created by Alejandro Quiroz on 3/20/11.
 *
 */

#include "OpMat.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
const double pi = 3.141592653589793238462643383279502884197;

/*============================================
 ========= Random Normal Multivariada =========
 ==============================================*/
// Aqui Y es de tamanio (d x n).
void Crmnorm(int n,int d, double *medias, double **Sg, double **TY)
{
	int i,j;
	double *AuxMu;
	double **AuxMu2, **Z, **Mu;
	double **TSg,**TZ,**Y;
	
	AuxMu=array1(n*d);
	
	TSg=array2(d,d);
	Z=array2(n,d);
	TZ=array2(d,n);
	Y=array2(d,n);
	
	
	GetRNGstate();	
	for (i=0;i<n*d;i++)
		AuxMu[i]=rnorm(0,1);
	PutRNGstate();
	
	AuxMu2=array2srce(AuxMu, n, d);
	
	MatChol(Sg,d);
	MatTrans(Sg, d, d,TSg);
	
	MatMult(AuxMu2, TSg,n, d, d,Z);
	MatTrans(Z, n, d,TZ);
	
	Mu=array2(d,n);
	for(i=0;i<n;i++)
	{
		for(j=0;j<d;j++)
		{
			Mu[j][i]=medias[j];
		}
	}
	
	MatSum(Mu,TZ, d, n,Y);
	MatTrans(Y, d, n,TY);
	
}
/*==============================================
 ========= Random Wishart Multivariada =========
 ===============================================*/

void Crwish(int v,int p, double **Sg, double **Aux3)
{
	int i,j,ii;
	double **Z, **Aux, **Aux2;
	double **TSg;
	
	TSg=array2(p,p);
	Z=array2(p,p);
	
	Aux=array2(p,p);
	Aux2=array2(p,p);
	
	for(i=0;i<p;i++)
		for(j=0;j<p;j++)
			Z[i][j]=0;
	
	GetRNGstate();
	for(i=0;i<p;i++)
	{
		ii=v-i;
		Z[i][i]=sqrt(rchisq(ii));
	}
	PutRNGstate();
	
	GetRNGstate();
	for(i=1;i<p;i++)
	{
		for(j=0;j<i;j++)
		{
			Z[j][i]=rnorm(0,1);
		}
	}
	PutRNGstate();
	MatChol(Sg,p);
	MatTrans(Sg, p, p, TSg);
	MatMult(Z, TSg, p, p, p,Aux);
	
	MatTrans(Aux, p, p,Aux2);
	
	MatMult(Aux2,Aux, p, p, p,Aux3);
	
}

/*============================================
 ========= Random Normal Multivariada =========
 ==============================================*/


double Cdmnorm(int *dd, double *xx, double *mean, double **Sg1, double **Sg2)
{
	int d,i,j;
	int *sQP;
	double suma, logDet,logPDF;
	double *X, *cQP, *dQP;
	double **Xvec, **ISg, **Q;
	
	d=*dd;
	sQP=arrayInt1(d);
	X=array1(d);
	cQP=array1(d);
	dQP=array1(d);
	
	ISg=array2(d,d);
	Q=array2(d,1);
	
	for(i=0;i<d;i++)
		X[i]=xx[i]-mean[i];
	Xvec=array2srce(X, d, 1);
	
	MatInv(Sg1,dd,ISg);	
	MatMult(ISg, Xvec, d, d, 1,Q);
	
	suma=0;
	for(j=0;j<d;j++)
		suma=suma+Xvec[j][0]*Q[j][0];
	
	MatQRdcmp(Sg2, d, cQP, dQP, sQP);
	logDet=0;
	for(i=0;i<d;i++)
		logDet=logDet+log(fabs(dQP[i]));
	
	logPDF=(suma+d*log(2*pi)+logDet)/(-2);
	
	return(logPDF);
}

