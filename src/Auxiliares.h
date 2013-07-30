//
//  Auxiliares.h
//  GibbsCS
//
//  Created by Alejandro Quiroz on 1/9/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <math.h>


//=============== Para sigma

double SgmYAuxCS(int noCol, double GrpSz, double *ymu, double param, double l0, double v0)
{
    int i;
    double M,aux1,scala,rate,forma, resultado;
    
    M=0;
    for(i=0;i<noCol;i++)
    {
        aux1=(ymu[i]-param)*(ymu[i]-param);
        M=M+aux1;
    }
    
    scala=(GrpSz*M/2)+l0;
    rate=1/scala;
    forma=(noCol/2)+v0;
    
    GetRNGstate();
    resultado=rgamma(forma,rate);
    PutRNGstate();
    
    resultado=1/resultado;
    return resultado;
}

//=============== Para Beta

double BetaAuxCS(int noColA, int noColB, double GrpSz, double *ymuA, double *ymuB, double param, double sgmA, double sgmB)
{
    int i;
    double Ma,Mb,M2,Mu,aux1,CAest,CBest,Cest,resultado;
    
    Ma=0;
    for(i=0;i<noColA;i++)
    {
        aux1=ymuA[i];
        Ma=Ma+aux1;
    }
    Ma=Ma/sgmA;
    
    Mb=0;
    for(i=0;i<noColB;i++)
    {
        aux1=ymuB[i]-param;
        Mb=Mb+aux1;
    }
    Mb=Mb/sgmB;
    
    M2=Ma+Mb;
    CAest=noColA*GrpSz/sgmA;
    CBest=noColB*GrpSz/sgmB;
    Cest=CAest+CBest;
    Cest=1/Cest;
    Mu=GrpSz*Cest*M2;
    Cest=sqrt(Cest);
    
    GetRNGstate();
    resultado=rnorm(Mu,Cest);
    PutRNGstate();
    
    return resultado;
}

//=============== Para Alfa

double AlfaAuxCS(int noColB, double GrpSz, double *ymuB, double param, double sgmB, double sgm, double rho, double mm)
{
    int i;
    double CSiDiff,CNoDiff,Cest,Mu,Mu2,aux1,aux2,aux3,aux4,auxSi,auxNo,suma,zuma,cons,rznLog,rzn,prob,alea,resultado;
    
    CSiDiff=(sgmB/GrpSz)+sgm;
    auxSi=sqrt(CSiDiff);
    CNoDiff=(sgmB/GrpSz);
    auxNo=sqrt(CNoDiff);
    
    suma=0;
    zuma=0;
    cons=0;
    Mu=0;
    for(i=0;i<noColB;i++)
    {
        aux1=ymuB[i]-param;
        aux2=dnorm(aux1,0,auxSi,1);
        aux3=dnorm(aux1,0,auxNo,1);
        
        suma=suma+aux2;
        zuma=zuma+aux3;
        Mu=Mu+aux1;
    }
    aux4=(rho*mm)/(1-rho*mm);
    cons=log(aux4);
    rznLog=suma-zuma+cons;

    if(rznLog>600){
        prob=1;
    }
    else if(rznLog<(-600)){
        prob=0;
    }
    else{
        rzn=exp(rznLog);
        prob=rzn/(1+rzn);
    }
    
    Cest=(noColB*GrpSz/sgmB)+(1/sgm);
    Cest=1/Cest;
    Mu2=(Mu*GrpSz*Cest/sgmB);
    
    GetRNGstate();
    alea=runif(0,1);
    PutRNGstate();
    
    resultado=0;
    if(alea<prob){
        Cest=sqrt(Cest);
        GetRNGstate();
        resultado=rnorm(Mu2,Cest);
        PutRNGstate();
    }
    return resultado;
}


