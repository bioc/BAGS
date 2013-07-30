//
//  Funciones.h
//  GibbsCS
//
//  Created by Alejandro Quiroz on 1/9/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Auxiliares.h"
#include "libmat.h"

//=============== Para sigma

double SgmYFunCS(double *noCol,double GrpSzs,double *YMu,double *Y2Mu,double B0,double A1,double A2,double A3,double A4,double L0, double V0)
{
    double forma;
    double escala1,escala2,escala3,escala,rate,resI;
    
    forma=((noCol[0]+noCol[1]+noCol[2]+noCol[3]+noCol[4])/2)+L0;
    
    escala1=(Y2Mu[0]+Y2Mu[1]+Y2Mu[2]+Y2Mu[3]+Y2Mu[4]);
    escala2=2*(noCol[0]*YMu[0]*B0+noCol[1]*YMu[1]*(B0+A1)+noCol[2]*YMu[2]*(B0+A2)+noCol[3]*YMu[3]*(B0+A3)+noCol[4]*YMu[4]*(B0+A4));
    escala3=(noCol[0]*B0*B0)+(noCol[1]*(B0+A1)*(B0+A1))+(noCol[2]*(B0+A2)*(B0+A2))+(noCol[3]*(B0+A3)*(B0+A3))+(noCol[4]*(B0+A4)*(B0+A4));
    escala=0.5*GrpSzs*(escala1-escala2+escala3)+V0;
    rate=1/escala;
    
    GetRNGstate();
    resI=rgamma(forma,rate);
    PutRNGstate();
    
    resI=1/resI;
    return resI;
}

double SgmYFunCSDiff(double noCol,double GrpSzs,double YMu,double Param,double L0, double V0)
{
    double forma;
    double escala1,escala,rate,resI;
    
    forma=(1/2);
    forma=forma+L0;
    
    escala1=(YMu-Param)*(YMu-Param)*noCol*GrpSzs*0.5;
    escala=escala1+V0;
    
    rate=1/escala;
    
    GetRNGstate();
    resI=rgamma(forma,rate);
    PutRNGstate();
    
    resI=1/resI;
    return resI;
}

double BetaFunCS(double *noCol,double *YMu,double GrpSzs,double A1,double A2,double A3,double A4,double SgmY)
{
    double media1,media,varianza1,varianza,resultado;
    
    media1=(noCol[0]*YMu[0])+noCol[1]*(YMu[1]-A1)+noCol[2]*(YMu[2]-A2)+noCol[3]*(YMu[3]-A3)+noCol[4]*(YMu[4]-A4);
    varianza1=noCol[0]+noCol[1]+noCol[2]+noCol[3]+noCol[4];
    
    media=media1/varianza1;
    varianza=SgmY/(varianza1*GrpSzs);
    varianza=sqrt(varianza);
    
    GetRNGstate();
    resultado=rnorm(media,varianza);
    PutRNGstate();
    
    return resultado;
}

double BetaFunCSDiff5(double *noCol,double *YMu,double GrpSzs,double A1,double A2,double A3,double A4,double SgmY0,double SgmY1,double SgmY2,double SgmY3,double SgmY4)
{
    double media1,media,varianza1,varianza,resultado,varaux0,varaux1,varaux2,varaux3,varaux4;
    
    varaux0=SgmY0;
    varaux1=SgmY1;
    varaux2=SgmY2;
    varaux3=SgmY3;
    varaux4=SgmY4;
    
    media1=(noCol[0]*GrpSzs*YMu[0]/varaux0)+(noCol[1]*GrpSzs*(YMu[1]-A1)/varaux1)+(noCol[2]*GrpSzs*(YMu[2]-A2)/varaux2)+(noCol[3]*GrpSzs*(YMu[3]-A3)/varaux3)+(noCol[4]*GrpSzs*(YMu[4]-A4)/varaux4);
    varianza1=(noCol[1]*GrpSzs/varaux1)+(noCol[1]*GrpSzs/varaux1)+(noCol[2]*GrpSzs/varaux2)+(noCol[3]*GrpSzs/varaux3)+(noCol[4]*GrpSzs/varaux4);
    
    media=(media1/varianza1);
    varianza=1/varianza1;
    varianza=sqrt(varianza);
    
    GetRNGstate();
    resultado=rnorm(media,varianza);
    PutRNGstate();
    
    return resultado;
}

double BetaFunCSDiff4(double *noCol,double *YMu,double GrpSzs,double A1,double A2,double A3,double SgmY0,double SgmY1,double SgmY2,double SgmY3)
{
    double media1,media,varianza1,varianza,resultado,varaux0,varaux1,varaux2,varaux3;
    
    varaux0=SgmY0;
    varaux1=SgmY1;
    varaux2=SgmY2;
    varaux3=SgmY3;
    
    media1=(noCol[0]*GrpSzs*YMu[0]/varaux0)+(noCol[1]*GrpSzs*(YMu[1]-A1)/varaux1)+(noCol[2]*GrpSzs*(YMu[2]-A2)/varaux2)+(noCol[3]*GrpSzs*(YMu[3]-A3)/varaux3);
    varianza1=(noCol[1]*GrpSzs/varaux1)+(noCol[1]*GrpSzs/varaux1)+(noCol[2]*GrpSzs/varaux2)+(noCol[3]*GrpSzs/varaux3);
    
    media=(media1/varianza1);
    varianza=1/varianza1;
    varianza=sqrt(varianza);
    
    GetRNGstate();
    resultado=rnorm(media,varianza);
    PutRNGstate();
    
    return resultado;
}

double BetaFunCSDiff3(double *noCol,double *YMu,double GrpSzs,double A1,double A2,double SgmY0,double SgmY1,double SgmY2)
{
    double media1,media,varianza1,varianza,resultado,varaux0,varaux1,varaux2;
    
    varaux0=SgmY1;
    varaux1=SgmY1;
    varaux2=SgmY2;
    
    media1=(noCol[0]*GrpSzs*YMu[0]/varaux0)+(noCol[1]*GrpSzs*(YMu[1]-A1)/varaux1)+(noCol[2]*GrpSzs*(YMu[2]-A2)/varaux2);
    varianza1=(noCol[1]*GrpSzs/varaux1)+(noCol[1]*GrpSzs/varaux1)+(noCol[2]*GrpSzs/varaux2);
    
    media=(media1/varianza1);
    varianza=1/varianza1;
    varianza=sqrt(varianza);
    
    GetRNGstate();
    resultado=rnorm(media,varianza);
    PutRNGstate();
    
    return resultado;
}

double BetaFunCSDiff2(double *noCol,double *YMu,double GrpSzs,double A1,double SgmY0,double SgmY1)
{
    double media1,media,varianza1,varianza,resultado,varaux0,varaux1,YMu1;
    
    varaux0=SgmY0;
    varaux1=SgmY1;
    
    media1=(noCol[0]*GrpSzs*YMu[0]/varaux0)+(noCol[1]*GrpSzs*(YMu[1]-A1)/varaux1);
    varianza1=(noCol[0]*GrpSzs/varaux0)+(noCol[1]*GrpSzs/varaux1);
    
    media=(media1/varianza1);
    varianza=1/varianza1;
    varianza=sqrt(varianza);
    
    GetRNGstate();
    resultado=rnorm(media,varianza);
    PutRNGstate();
    
    return resultado;
}

double AlfaFunCS(double noCol,double YMu,double GrpSzs,double B0,double SgmY,double SgmA,double rho,double mm)
{
    double CSiDiff,CNoDiff,Cest,Mu,aux1,aux4,suma,zuma,cons,rznLog,rzn,prob,alea,resultado,C,A,varianza,MediaDiff;
    
    CSiDiff=(SgmY/(noCol*GrpSzs))+SgmA;
    CSiDiff=sqrt(CSiDiff);
    CNoDiff=(SgmY/(noCol*GrpSzs));
    CNoDiff=sqrt(CNoDiff);
    
    aux1=YMu-B0;
    suma=dnorm(aux1,0,CSiDiff,1);
    zuma=dnorm(aux1,0,CNoDiff,1);
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
    
    Cest=(GrpSzs*noCol/SgmY)+(1/SgmA);
    varianza=(1/Cest);
    Mu=GrpSzs*noCol*(YMu-B0)/SgmY;
    Mu=Mu/Cest;
    
    GetRNGstate();
    alea=runif(0,1);
    PutRNGstate();
    
    resultado=0;
    if(alea<prob){
        varianza=sqrt(varianza);
        GetRNGstate();
        resultado=rnorm(Mu,varianza);
        PutRNGstate();
    }
    return resultado;
}

