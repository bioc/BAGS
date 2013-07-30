//
//  Created by Alejandro Quiroz on 3/4/13.
//
// 4
// 5 Created by Alejandro Quiroz on 1/9/12.
// 7

#include "FuncionesFast-Final.h"

void Gibbs5(int *noRow,double *noCol,int *iter,double *GrpSzs,double *Yu,double *L0,double *V0,double *L0A,double *V0A,double *MM,double *AAPi,double *ApriDiffExp,double *resultado1,double *resultado2,double *resultado3,double *resultado4)
{
    int i,j,nR,nI,contador;
    double *nC;
    double *YMu;
    double Mm,AaPi,Ade,Ss,Rr,c0,prob,alea,forma1,forma2,shape1,shape2;
    double **YMuMat,*forma,*scala,*formaA,*scalaA;
    
    double ResRho1,ResRho2,ResRho3,ResRho4,ResSgmA1,ResSgmA2,ResSgmA3,ResSgmA4,ResSgmA;
    double *ResA1,*ResA2,*ResA3,*ResA4,*ResB,*ResSgmY0,*ResSgmY1,*ResSgmY2,*ResSgmY3,*ResSgmY4,*ResSgmY;
    
    double suma_ro1,suma_ro2,suma_ro3,suma_ro4,suma_sa1,suma_sa2,suma_sa3,suma_sa4,suma_no1,suma_no2,suma_no3,suma_no4,parametro;
    
    //
    double formaAux;
    double escala1,escala2,escala3,escala,rate,resI;
    //
    
    nR=*noRow;
    nI=*iter;
    Mm=*MM;
    AaPi=*AAPi;
    Ade=*ApriDiffExp;
    
    nC=array1srce(noCol,5);
    YMuMat=array2srce(Yu,nR,5);
    YMu=array1(5);
    
    forma=array1srce(L0,5);
    scala=array1srce(V0,5);
    formaA=array1srce(L0A,4);
    scalaA=array1srce(V0A,4);
    
    ResSgmA=0.1;
    ResSgmA1=0.1;
    ResSgmA2=0.1;
    ResSgmA3=0.1;
    ResSgmA4=0.1;
    ResRho1=0.4;
    ResRho2=0.4;
    ResRho3=0.4;
    ResRho4=0.4;
    
    ResSgmY=array1(nR);
    ResSgmY0=array1(nR);
    ResSgmY1=array1(nR);
    ResSgmY2=array1(nR);
    ResSgmY3=array1(nR);
    ResSgmY4=array1(nR);
    ResB=array1(nR);
    ResA1=array1(nR);
    ResA2=array1(nR);
    ResA3=array1(nR);
    ResA4=array1(nR);
    
    Ss=nR;
    Rr=Ade/nR;
    for(j=0;j<nR;j++)
    {
        ResSgmY[j]=0;
        ResSgmY0[j]=0;
        ResSgmY1[j]=0;
        ResSgmY2[j]=0;
        ResSgmY3[j]=0;
        ResSgmY4[j]=0;
        
        ResB[j]=6;
        
        ResA1[j]=0.1;
        ResA2[j]=0.2;
        ResA3[j]=0.3;
        ResA4[j]=0.4;
    }
    contador=0;
    for(i=0;i<nI;i++)
    {
        suma_ro1=0;
        suma_ro2=0;
        suma_ro3=0;
        suma_ro4=0;
        
        suma_sa1=0;
        suma_sa2=0;
        suma_sa3=0;
        suma_sa4=0;
        
        suma_no1=0;
        suma_no2=0;
        suma_no3=0;
        suma_no4=0;
        for(j=0;j<nR;j++)
        {
            YMu[0]=YMuMat[j][0];
            YMu[1]=YMuMat[j][1];
            YMu[2]=YMuMat[j][2];
            YMu[3]=YMuMat[j][3];
            YMu[4]=YMuMat[j][4];
               
            ResSgmY0[j]=SgmYFunCSDiff(nC[0],GrpSzs[j],YMu[0],ResB[j],forma[0],scala[0]);
            
            parametro=ResB[j]+ResA1[j];
            ResSgmY1[j]=SgmYFunCSDiff(nC[1],GrpSzs[j],YMu[1],parametro,forma[1],scala[1]);
            
            parametro=ResB[j]+ResA2[j];
            ResSgmY2[j]=SgmYFunCSDiff(nC[2],GrpSzs[j],YMu[2],parametro,forma[2],scala[2]);
            
            parametro=ResB[j]+ResA3[j];
            ResSgmY3[j]=SgmYFunCSDiff(nC[3],GrpSzs[j],YMu[3],parametro,forma[3],scala[3]);
            
            parametro=ResB[j]+ResA4[j];
            ResSgmY4[j]=SgmYFunCSDiff(nC[4],GrpSzs[j],YMu[4],parametro,forma[4],scala[4]);
            
            ResB[j]=BetaFunCSDiff5(nC,YMu,GrpSzs[j],ResA1[j],ResA2[j],ResA3[j],ResA4[j],ResSgmY0[j],ResSgmY1[j],ResSgmY2[j],ResSgmY3[j],ResSgmY4[j]);
            ResA1[j]=AlfaFunCS(nC[1],YMu[1],GrpSzs[j],ResB[j],ResSgmY1[j],ResSgmA1,ResRho1,Mm);
            ResA2[j]=AlfaFunCS(nC[2],YMu[2],GrpSzs[j],ResB[j],ResSgmY2[j],ResSgmA2,ResRho2,Mm);
            ResA3[j]=AlfaFunCS(nC[3],YMu[3],GrpSzs[j],ResB[j],ResSgmY3[j],ResSgmA3,ResRho3,Mm);
            ResA4[j]=AlfaFunCS(nC[4],YMu[4],GrpSzs[j],ResB[j],ResSgmY4[j],ResSgmA4,ResRho4,Mm);
            
            //========== Parametro PI 1 ===========
            if(ResA1[j]==0)
            {
                c0=(1-Mm)*ResRho1/(1-ResRho1);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro1=suma_ro1+1;
                
                suma_no1=suma_no1+1;
            }
            else{
                suma_ro1=suma_ro1+1;
            }
            //========== Parametro PI 2 ===========
            if(ResA2[j]==0)
            {
                c0=(1-Mm)*ResRho2/(1-ResRho2);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro2=suma_ro2+1;
                
                suma_no2=suma_no2+1;
            }
            else{
                suma_ro2=suma_ro2+1;
            }
            //========== Parametro PI 3 ===========
            if(ResA3[j]==0)
            {
                c0=(1-Mm)*ResRho3/(1-ResRho3);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro3=suma_ro3+1;
                
                suma_no3=suma_no3+1;
            }
            else{
                suma_ro3=suma_ro3+1;
            }
            //========== Parametro PI 4 ===========
            if(ResA4[j]==0)
            {
                c0=(1-Mm)*ResRho4/(1-ResRho4);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro4=suma_ro4+1;
                
                suma_no4=suma_no4+1;
            }
            else{
                suma_ro4=suma_ro4+1;
            }
            //========== Parametro Sigma Alfa 1 ==========
            suma_sa1=suma_sa1+(ResA1[j]*ResA1[j]);
            //========== Parametro Sigma Alfa 2 ==========
            suma_sa2=suma_sa2+(ResA2[j]*ResA2[j]);
            //========== Parametro Sigma Alfa 3 ==========
            suma_sa3=suma_sa3+(ResA3[j]*ResA3[j]);
            //========== Parametro Sigma Alfa 4 ==========
            suma_sa4=suma_sa4+(ResA4[j]*ResA4[j]);
            
            resultado1[contador]=ResA1[j];
            resultado2[contador]=ResA2[j];
            resultado3[contador]=ResA3[j];
            resultado4[contador]=ResA4[j];
            contador=contador+1;
            
            //Rprintf("Calculo del Grupo:\t");
            //Rprintf("%d\n",j);
        }
        //========== Parametro Sigma Alfa 1 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no1)/2;
        forma1=forma1+formaA[0];
        //Rprintf("Forma1\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa1/2);
        forma2=forma2+scalaA[0];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA1=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA1=1/ResSgmA1;
        //========== Parametro Sigma Alfa 2 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no2)/2;
        forma1=forma1+formaA[1];
        //Rprintf("Forma2\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa2/2);
        forma2=forma2+scalaA[1];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA2=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA2=1/ResSgmA2;
        //========== Parametro Sigma Alfa 3 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no3)/2;
        forma1=forma1+formaA[2];
        //Rprintf("Forma3\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa3/2);
        forma2=forma2+scalaA[2];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA3=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA3=1/ResSgmA3;
        //========== Parametro Sigma Alfa 4 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no4)/2;
        forma1=forma1+formaA[3];
        //Rprintf("Forma4\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa4/2);
        forma2=forma2+scalaA[3];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA4=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA4=1/ResSgmA4;
        
        /*forma1=((suma_no1+suma_no2+suma_no3+suma_no4)/2)+formaA[0];
        forma2=((suma_sa1+suma_sa2+suma_sa3+suma_sa4)/2)+scalaA[0];
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA=1/ResSgmA;*/
        
        //========== Parametro Rho 1 ==========
        shape1=Ss*Rr+suma_ro1;
        shape2=Ss*(1-Rr)+(nR-suma_ro1);
        GetRNGstate();
        ResRho1=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 2 ==========
        shape1=Ss*Rr+suma_ro2;
        shape2=Ss*(1-Rr)+(nR-suma_ro2);
        GetRNGstate();
        ResRho2=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 3 ==========
        shape1=Ss*Rr+suma_ro3;
        shape2=Ss*(1-Rr)+(nR-suma_ro3);
        GetRNGstate();
        ResRho3=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 4 ==========
        shape1=Ss*Rr+suma_ro4;
        shape2=Ss*(1-Rr)+(nR-suma_ro4);
        GetRNGstate();
        ResRho4=rbeta(shape1,shape2);
        PutRNGstate();
        
        //Rprintf("FIN DE LA ITERACION:\t");
        //Rprintf("%d\n",i);
    }
    
}


void Gibbs4(int *noRow,double *noCol,int *iter,double *GrpSzs,double *Yu,double *L0,double *V0,double *L0A,double *V0A,double *MM,double *AAPi,double *ApriDiffExp,double *resultado1,double *resultado2,double *resultado3)
{
    int i,j,nR,nI,contador;
    double *nC;
    double *YMu;
    double Mm,AaPi,Ade,Ss,Rr,c0,prob,alea,forma1,forma2,shape1,shape2;
    double **YMuMat,*forma,*scala,*formaA,*scalaA;
    
    double ResRho1,ResRho2,ResRho3,ResSgmA1,ResSgmA2,ResSgmA3,ResSgmA;
    double *ResA1,*ResA2,*ResA3,*ResB,*ResSgmY0,*ResSgmY1,*ResSgmY2,*ResSgmY3,*ResSgmY;
    
    double suma_ro1,suma_ro2,suma_ro3,suma_sa1,suma_sa2,suma_sa3,suma_no1,suma_no2,suma_no3,parametro;
    
    //
    double formaAux;
    double escala1,escala2,escala3,escala,rate,resI;    // Aqui tengo duda
    //
    
    nR=*noRow;
    nI=*iter;
    Mm=*MM;
    AaPi=*AAPi;
    Ade=*ApriDiffExp;
    
    nC=array1srce(noCol,4);
    YMuMat=array2srce(Yu,nR,4);
    YMu=array1(4);

    forma=array1srce(L0,4);
    scala=array1srce(V0,4);
    formaA=array1srce(L0A,3);
    scalaA=array1srce(V0A,3);
    
    ResSgmA=0.1;
    ResSgmA1=0.1;
    ResSgmA2=0.1;
    ResSgmA3=0.1;
    ResRho1=0.4;
    ResRho2=0.4;
    ResRho3=0.4;
    
    ResSgmY=array1(nR);
    ResSgmY0=array1(nR);
    ResSgmY1=array1(nR);
    ResSgmY2=array1(nR);
    ResSgmY3=array1(nR);
    ResB=array1(nR);
    ResA1=array1(nR);
    ResA2=array1(nR);
    ResA3=array1(nR);
    
    Ss=nR;
    Rr=Ade/nR;
    for(j=0;j<nR;j++)
    {
        ResSgmY[j]=0;
        ResSgmY0[j]=0;
        ResSgmY1[j]=0;
        ResSgmY2[j]=0;
        ResSgmY3[j]=0;
        
        ResB[j]=6;
        
        ResA1[j]=0.1;
        ResA2[j]=0.2;
        ResA3[j]=0.3;
    }
    contador=0;
    for(i=0;i<nI;i++)
    {
        suma_ro1=0;
        suma_ro2=0;
        suma_ro3=0;
        
        suma_sa1=0;
        suma_sa2=0;
        suma_sa3=0;
        
        suma_no1=0;
        suma_no2=0;
        suma_no3=0;
        for(j=0;j<nR;j++)
        {
            YMu[0]=YMuMat[j][0];
            YMu[1]=YMuMat[j][1];
            YMu[2]=YMuMat[j][2];
            YMu[3]=YMuMat[j][3];
            
            ResSgmY0[j]=SgmYFunCSDiff(nC[0],GrpSzs[j],YMu[0],ResB[j],forma[0],scala[0]);
            
            parametro=ResB[j]+ResA1[j];
            ResSgmY1[j]=SgmYFunCSDiff(nC[1],GrpSzs[j],YMu[1],parametro,forma[1],scala[1]);
            
            parametro=ResB[j]+ResA2[j];
            ResSgmY2[j]=SgmYFunCSDiff(nC[2],GrpSzs[j],YMu[2],parametro,forma[2],scala[2]);
            
            parametro=ResB[j]+ResA3[j];
            ResSgmY3[j]=SgmYFunCSDiff(nC[3],GrpSzs[j],YMu[3],parametro,forma[3],scala[3]);
            
            ResB[j]=BetaFunCSDiff4(nC,YMu,GrpSzs[j],ResA1[j],ResA2[j],ResA3[j],ResSgmY0[j],ResSgmY1[j],ResSgmY2[j],ResSgmY3[j]);
            ResA1[j]=AlfaFunCS(nC[1],YMu[1],GrpSzs[j],ResB[j],ResSgmY1[j],ResSgmA1,ResRho1,Mm);
            ResA2[j]=AlfaFunCS(nC[2],YMu[2],GrpSzs[j],ResB[j],ResSgmY2[j],ResSgmA2,ResRho2,Mm);
            ResA3[j]=AlfaFunCS(nC[3],YMu[3],GrpSzs[j],ResB[j],ResSgmY3[j],ResSgmA3,ResRho3,Mm);
            
            //========== Parametro PI 1 ===========
            if(ResA1[j]==0)
            {
                c0=(1-Mm)*ResRho1/(1-ResRho1);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro1=suma_ro1+1;
                
                suma_no1=suma_no1+1;
            }
            else{
                suma_ro1=suma_ro1+1;
            }
            //========== Parametro PI 2 ===========
            if(ResA2[j]==0)
            {
                c0=(1-Mm)*ResRho2/(1-ResRho2);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro2=suma_ro2+1;
                
                suma_no2=suma_no2+1;
            }
            else{
                suma_ro2=suma_ro2+1;
            }
            //========== Parametro PI 3 ===========
            if(ResA3[j]==0)
            {
                c0=(1-Mm)*ResRho3/(1-ResRho3);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro3=suma_ro3+1;
                
                suma_no3=suma_no3+1;
            }
            else{
                suma_ro3=suma_ro3+1;
            }
            
            //========== Parametro Sigma Alfa 1 ==========
            suma_sa1=suma_sa1+(ResA1[j]*ResA1[j]);
            //========== Parametro Sigma Alfa 2 ==========
            suma_sa2=suma_sa2+(ResA2[j]*ResA2[j]);
            //========== Parametro Sigma Alfa 3 ==========
            suma_sa3=suma_sa3+(ResA3[j]*ResA3[j]);
            
            resultado1[contador]=ResA1[j];
            resultado2[contador]=ResA2[j];
            resultado3[contador]=ResA3[j];
            contador=contador+1;
            
            //Rprintf("Calculo del Grupo:\t");
            //Rprintf("%d\n",j);
        }
        //========== Parametro Sigma Alfa 1 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no1)/2;
        forma1=forma1+formaA[0];
        //Rprintf("Forma1\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa1/2);
        forma2=forma2+scalaA[0];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA1=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA1=1/ResSgmA1;
        //========== Parametro Sigma Alfa 2 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no2)/2;
        forma1=forma1+formaA[1];
        //Rprintf("Forma2\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa2/2);
        forma2=forma2+scalaA[1];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA2=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA2=1/ResSgmA2;
        //========== Parametro Sigma Alfa 3 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no3)/2;
        forma1=forma1+formaA[2];
        //Rprintf("Forma3\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa3/2);
        forma2=forma2+scalaA[2];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA3=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA3=1/ResSgmA3;
        
        //========== Parametro Rho 1 ==========
        shape1=Ss*Rr+suma_ro1;
        shape2=Ss*(1-Rr)+(nR-suma_ro1);
        GetRNGstate();
        ResRho1=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 2 ==========
        shape1=Ss*Rr+suma_ro2;
        shape2=Ss*(1-Rr)+(nR-suma_ro2);
        GetRNGstate();
        ResRho2=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 3 ==========
        shape1=Ss*Rr+suma_ro3;
        shape2=Ss*(1-Rr)+(nR-suma_ro3);
        GetRNGstate();
        ResRho3=rbeta(shape1,shape2);
        PutRNGstate();
        //Rprintf("FIN DE LA ITERACION:\t");
        //Rprintf("%d\n",i);
    }
    
    
}

void Gibbs3(int *noRow,double *noCol,int *iter,double *GrpSzs,double *Yu,double *L0,double *V0,double *L0A,double *V0A,double *MM,double *AAPi,double *ApriDiffExp,double *resultado1,double *resultado2)
{
    int i,j,nR,nI,contador;
    double *nC;
    double *YMu;
    double Mm,AaPi,Ade,Ss,Rr,c0,prob,alea,forma1,forma2,shape1,shape2;
    double **YMuMat,*forma,*scala,*formaA,*scalaA;
    
    double ResRho1,ResRho2,ResSgmA1,ResSgmA2,ResSgmA;
    double *ResA1,*ResA2,*ResB,*ResSgmY0,*ResSgmY1,*ResSgmY2,*ResSgmY;
    
    double suma_ro1,suma_ro2,suma_sa1,suma_sa2,suma_no1,suma_no2,parametro;
    
    //
    double formaAux;
    double escala1,escala2,escala3,escala,rate,resI;    // Aqui tengo duda
    //
    
    nR=*noRow;
    nI=*iter;
    Mm=*MM;
    AaPi=*AAPi;
    Ade=*ApriDiffExp;
    
    nC=array1srce(noCol,3);
    YMuMat=array2srce(Yu,nR,3);
    YMu=array1(3);

    forma=array1srce(L0,3);
    scala=array1srce(V0,3);
    formaA=array1srce(L0A,2);
    scalaA=array1srce(V0A,2);
    
    ResSgmA=0.1;
    ResSgmA1=0.1;
    ResSgmA2=0.1;
    ResRho1=0.4;
    ResRho2=0.4;
    
    ResSgmY=array1(nR);
    ResSgmY0=array1(nR);
    ResSgmY1=array1(nR);
    ResSgmY2=array1(nR);
    ResB=array1(nR);
    ResA1=array1(nR);
    ResA2=array1(nR);
    
    Ss=nR;
    Rr=Ade/nR;
    for(j=0;j<nR;j++)
    {
        ResSgmY[j]=0;
        ResSgmY0[j]=0;
        ResSgmY1[j]=0;
        ResSgmY2[j]=0;
        
        ResB[j]=6;
        
        ResA1[j]=0.1;
        ResA2[j]=0.2;
    }
    contador=0;
    for(i=0;i<nI;i++)
    {
        suma_ro1=0;
        suma_ro2=0;
        
        suma_sa1=0;
        suma_sa2=0;
        
        suma_no1=0;
        suma_no2=0;
        for(j=0;j<nR;j++)
        {
            YMu[0]=YMuMat[j][0];
            YMu[1]=YMuMat[j][1];
            YMu[2]=YMuMat[j][2];
            
            ResSgmY0[j]=SgmYFunCSDiff(nC[0],GrpSzs[j],YMu[0],ResB[j],forma[0],scala[0]);
            
            parametro=ResB[j]+ResA1[j];
            ResSgmY1[j]=SgmYFunCSDiff(nC[1],GrpSzs[j],YMu[1],parametro,forma[1],scala[1]);
            
            parametro=ResB[j]+ResA2[j];
            ResSgmY2[j]=SgmYFunCSDiff(nC[2],GrpSzs[j],YMu[2],parametro,forma[2],scala[2]);
            
            ResB[j]=BetaFunCSDiff3(nC,YMu,GrpSzs[j],ResA1[j],ResA2[j],ResSgmY0[j],ResSgmY1[j],ResSgmY2[j]);
            ResA1[j]=AlfaFunCS(nC[1],YMu[1],GrpSzs[j],ResB[j],ResSgmY1[j],ResSgmA1,ResRho1,Mm);
            ResA2[j]=AlfaFunCS(nC[2],YMu[2],GrpSzs[j],ResB[j],ResSgmY2[j],ResSgmA2,ResRho2,Mm);
            
            //========== Parametro PI 1 ===========
            if(ResA1[j]==0)
            {
                c0=(1-Mm)*ResRho1/(1-ResRho1);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro1=suma_ro1+1;
                
                suma_no1=suma_no1+1;
            }
            else{
                suma_ro1=suma_ro1+1;
            }
            //========== Parametro PI 2 ===========
            if(ResA2[j]==0)
            {
                c0=(1-Mm)*ResRho2/(1-ResRho2);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro2=suma_ro2+1;
                
                suma_no2=suma_no2+1;
            }
            else{
                suma_ro2=suma_ro2+1;
            }
            
            //========== Parametro Sigma Alfa 1 ==========
            suma_sa1=suma_sa1+(ResA1[j]*ResA1[j]);
            //========== Parametro Sigma Alfa 2 ==========
            suma_sa2=suma_sa2+(ResA2[j]*ResA2[j]);
            
            resultado1[contador]=ResA1[j];
            resultado2[contador]=ResA2[j];
            contador=contador+1;
            
            //Rprintf("Calculo del Grupo:\t");
            //Rprintf("%d\n",j);
        }
        //========== Parametro Sigma Alfa 1 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no1)/2;
        forma1=forma1+formaA[0];
        //Rprintf("Forma1\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa1/2);
        forma2=forma2+scalaA[0];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA1=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA1=1/ResSgmA1;
        //========== Parametro Sigma Alfa 2 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no2)/2;
        forma1=forma1+formaA[1];
        //Rprintf("Forma2\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa2/2);
        forma2=forma2+scalaA[1];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA2=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA2=1/ResSgmA2;
        
        //========== Parametro Rho 1 ==========
        shape1=Ss*Rr+suma_ro1;
        shape2=Ss*(1-Rr)+(nR-suma_ro1);
        GetRNGstate();
        ResRho1=rbeta(shape1,shape2);
        PutRNGstate();
        //========== Parametro Rho 2 ==========
        shape1=Ss*Rr+suma_ro2;
        shape2=Ss*(1-Rr)+(nR-suma_ro2);
        GetRNGstate();
        ResRho2=rbeta(shape1,shape2);
        PutRNGstate();
        //Rprintf("FIN DE LA ITERACION:\t");
        //Rprintf("%d\n",i);
    }
    
    
}

void Gibbs2(int *noRow,double *noCol,int *iter,double *GrpSzs,double *Yu,double *L0,double *V0,double *L0A,double *V0A,double *MM,double *AAPi,double *ApriDiffExp,double *resultado1)
{
    int i,j,nR,nI,contador;
    double *nC;
    double *YMu;
    double Mm,AaPi,Ade,Ss,Rr,c0,prob,alea,forma1,forma2,shape1,shape2;
    double **YMuMat,*forma,*scala,*formaA,*scalaA;
    
    double ResRho1,ResSgmA1,ResSgmA;
    double *ResA1,*ResB,*ResSgmY0,*ResSgmY1,*ResSgmY;
    
    double suma_ro1,suma_sa1,suma_no1,parametro;
    
    //
    double formaAux;
    double escala1,escala2,escala3,escala,rate,resI;    // Aqui tengo duda
    //
    
    nR=*noRow;
    nI=*iter;
    Mm=*MM;
    AaPi=*AAPi;
    Ade=*ApriDiffExp;
    
    nC=array1srce(noCol,2);
    YMuMat=array2srce(Yu,nR,2);
    YMu=array1(2);
    
    forma=array1srce(L0,2);
    scala=array1srce(V0,2);
    formaA=array1srce(L0A,1);
    scalaA=array1srce(V0A,1);
    
    ResSgmA=0.1;
    ResSgmA1=0.1;
    ResRho1=0.4;
    
    ResSgmY=array1(nR);
    ResSgmY0=array1(nR);
    ResSgmY1=array1(nR);
    ResB=array1(nR);
    ResA1=array1(nR);
    
    Ss=nR;
    Rr=Ade/nR;
    for(j=0;j<nR;j++)
    {
        ResSgmY[j]=0;
        ResSgmY0[j]=0;
        ResSgmY1[j]=0;
        
        ResB[j]=6;
        
        ResA1[j]=0.1;
    }
    contador=0;
    for(i=0;i<nI;i++)
    {
        suma_ro1=0;
        
        suma_sa1=0;
        
        suma_no1=0;
        for(j=0;j<nR;j++)
        {
            YMu[0]=YMuMat[j][0];
            YMu[1]=YMuMat[j][1];
            
            ResSgmY0[j]=SgmYFunCSDiff(nC[0],GrpSzs[j],YMu[0],ResB[j],forma[0],scala[0]);
            
            parametro=ResB[j]+ResA1[j];
            ResSgmY1[j]=SgmYFunCSDiff(nC[1],GrpSzs[j],YMu[1],parametro,forma[1],scala[1]);
            
            ResB[j]=BetaFunCSDiff2(nC,YMu,GrpSzs[j],ResA1[j],ResSgmY0[j],ResSgmY1[j]);
            ResA1[j]=AlfaFunCS(nC[1],YMu[1],GrpSzs[j],ResB[j],ResSgmY1[j],ResSgmA1,ResRho1,Mm);
            
            //========== Parametro PI 1 ===========
            if(ResA1[j]==0)
            {
                c0=(1-Mm)*ResRho1/(1-ResRho1);
                prob=c0/(1+c0);
                
                GetRNGstate();
                alea=runif(0,1);
                PutRNGstate();
                
                if(alea<prob)
                    suma_ro1=suma_ro1+1;
                
                suma_no1=suma_no1+1;
            }
            else{
                suma_ro1=suma_ro1+1;
            }
            
            //========== Parametro Sigma Alfa 1 ==========
            suma_sa1=suma_sa1+(ResA1[j]*ResA1[j]);
            
            resultado1[contador]=ResA1[j];
            contador=contador+1;
            
            //Rprintf("Calculo del Grupo:\t");
            //Rprintf("%d\n",j);
        }
        //========== Parametro Sigma Alfa 1 ==========
        //forma1=(suma_no1/2)+formaA[0];
        forma1=(nR-suma_no1)/2;
        forma1=forma1+formaA[0];
        //Rprintf("Forma1\t");
        //Rprintf("%f\n",forma1);
        
        //forma2=(suma_sa1/2)+scalaA[0];
        forma2=(suma_sa1/2);
        forma2=forma2+scalaA[0];
        //Rprintf("%f\n",forma2);
        
        forma2=1/forma2;
        GetRNGstate();
        ResSgmA1=rgamma(forma1,forma2);
        PutRNGstate();
        ResSgmA1=1/ResSgmA1;
        
        //========== Parametro Rho 1 ==========
        shape1=Ss*Rr+suma_ro1;
        shape2=Ss*(1-Rr)+(nR-suma_ro1);
        GetRNGstate();
        ResRho1=rbeta(shape1,shape2);
        PutRNGstate();
        //Rprintf("FIN DE LA ITERACION:\t");
        //Rprintf("%d\n",i);
    }
    
    
}
