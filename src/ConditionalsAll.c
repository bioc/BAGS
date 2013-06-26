//
//  ConditionalsAll.c
//
//  Created by Alejandro Quiroz on 10/5/11.
//

#include "ConditionalsAll.h"
#include "InvWishart.h"


R_CMethodDef cMethods[] = {
    {"conditionalSg", (DL_FUNC) &conditionalSg, 12},
    {"conditionalBeta", (DL_FUNC) &conditionalBeta, 10},
    {"conditionalNPAlfa0", (DL_FUNC) &conditionalNPAlfa0, 10},
    {"conditionalNPAlfa", (DL_FUNC) &conditionalNPAlfa, 13},
    {"conditionalALLAlfa", (DL_FUNC) &conditionalALLAlfa, 12},
    {NULL, NULL, 0}
};

//==================================================================

void R_init_GeneGroupAnalysis(DllInfo *info) {
    R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}


//==================================================================
//==========    Sigma_i    =========================================
//==================================================================

void conditionalSg(int *II,int *TT,int *DD,int *indices,double *TamGrp,double *YA,double *YB,double *betaI,double *alfaI,double *LL0,double *vv0,double *Resul)
{
	int i,j,k,h,a,contador;
	int I,T,D,IT;
	int *ind;
	double v0,vn,ki;
	double **TamGrupos,**Beta,**Alfa,**YmA,**YmB,**L0;
	double **M,**AuxM,**MuA,**MuB,**MA,**MB,**MC,**SubYmA,**SubYmB,**SubYmAi,**SubYmBi,**rSubYmAi,**rSubYmBi,**TSubYmAi,**TSubYmBi,**Ln;
	double **AuxLn,**RLn,**Resultado;
	
	/* ========================================
	 /*  Definicion de las variables de trabajo
	 ======================================== */
	
	// Ints
	I=*II;
	T=*TT;
	D=*DD;
	IT=I*T;
	
	//Doubles
	v0=*vv0;
	
	// Matrices con datos
	YmA=array2srce(YA,I,T*D); //----------------- Liberar memoria
	YmB=array2srce(YB,I,T*D); //----------------- Liberar memoria
	TamGrupos=array2srce(TamGrp, I, 1); //------- Liberar memoria
	Beta=array2srce(betaI, I, T); //------------- Liberar memoria
	Alfa=array2srce(alfaI, I, T); //------------- Liberar memoria
	L0=array2srce(LL0, T, T); //----------------- Liberar memoria
	
	
	// Matrices incializadas en 0's
	M=array2(T,T); //---------------------------- Liberar memoria
	AuxM=array2(T,T); //------------------------- Liberar memoria
	MuA=array2(T,1);  //------------------------- Liberar memoria
	MuB=array2(T,1);  //------------------------- Liberar memoria
	ind=arrayInt1(T); //------------------------- Liberar memoria
	MA=array2(T,T); //--------------------------- Liberar memoria
	MB=array2(T,T); //--------------------------- Liberar memoria
	MC=array2(T,T); //--------------------------- Liberar memoria
	
	SubYmA=array2(I,T);  //---------------------- Liberar memoria
	SubYmB=array2(I,T);  //---------------------- Liberar memoria
	SubYmAi=array2(T,1); //---------------------- Liberar memoria
	SubYmBi=array2(T,1); //---------------------- Liberar memoria
	rSubYmAi=array2(T,1); //--------------------- Liberar memoria
	rSubYmBi=array2(T,1); //--------------------- Liberar memoria
	TSubYmAi=array2(1,T); //--------------------- Liberar memoria
	TSubYmBi=array2(1,T); //--------------------- Liberar memoria
	
	Ln=array2(T,T); //--------------------------- Liberar memoria
	AuxLn=array2(T,T); //------------------------ Liberar memoria 
	RLn=array2(T,T); //-------------------------- Liberar memoria
	
	Resultado=array2(I*T,T); //------------------ Liberar memoria
	
	vn=(2*D)+v0;
	
	for(i=0;i<I;i++)
	{
		ki=TamGrupos[i][0];
		
		// Poniendo 0's la matriz M
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				M[j][k]=0;
			}
		}
		
		
		for(j=0;j<D;j++)
		{
			// Definicion del vector de indices y de las medias de los grupos
			for(k=0;k<T;k++)
			{
				ind[k]=indices[k]+j;
				MuA[k][0]=Beta[i][k];
				MuB[k][0]=Beta[i][k]+Alfa[i][k];
			}
			
			SubMat(YmA, I, ind, T,SubYmA);
			SubMat(YmB, I, ind, T,SubYmB);
			for(k=0;k<T;k++)
			{
				SubYmAi[k][0]=SubYmA[i][k];
				SubYmBi[k][0]=SubYmB[i][k];
			}
			MatRes(SubYmAi,MuA,T,1,rSubYmAi);
			MatRes(SubYmBi,MuB,T,1,rSubYmBi);
			
			MatTrans(rSubYmAi,T,1,TSubYmAi);
			MatTrans(rSubYmBi,T,1,TSubYmBi);
			
			MatMult(rSubYmAi,TSubYmAi,T,1,T,MA);
			MatMult(rSubYmBi,TSubYmBi,T,1,T,MB);
			MatSum(MA, MB, T, T,MC);
			MatSum(M, MC, T, T,AuxM);
			
			for(k=0;k<T;k++)
			{
				for(h=0;h<T;h++)
				{
					M[k][h]=AuxM[k][h];
					AuxM[k][h]=0;
				}
			}
		}
		
		for(k=0;k<T;k++)
		{
			for(h=0;h<T;h++)
			{
				M[k][h]=ki*M[k][h];
			}
		}
		
		MatSum(M, L0, T, T, Ln);
		MatInv(Ln,TT,AuxLn);
		Criwish(vn, TT, AuxLn,RLn);
		
		a=T*i;
		for(k=0;k<T;k++)
		{
			for(h=0;h<T;h++)
			{
				Resultado[a+h][k]=RLn[h][k];
			}
		}
	}
	
	contador=0;
	for(k=0;k<T;k++)
	{
		for(j=0;j<IT;j++)
		{
			Resul[contador]=Resultado[j][k];
			contador=contador+1;
		}
	}
	
}

//===============================================================
//==========    Beta    =========================================
//===============================================================

void conditionalBeta(int *II,int *TT,int *DD,int *indices,double *TamGrp,double *YA,double *YB,double *alfaI,double *Si,double *Resul)
{
	int i,k,j,a,contador;
	int I,T,D,IT;
	int *ind;
	double ki;
	double *AuxMuFinal3;
	double **TamGrupos,**Alfa,**YmA,**YmB,**SSi;
	double **Mu,**AuxMu,**AuxMu2,**MuA,**MuB,**SubYmA,**SubYmB,**SubYmAi,**SubYmBi,**rSubYmAi,**rSubYmBi;
	double **Si1,**InvSi,**Cest,**InvCest,**AuxMuFinal,**AuxMuFinal2;
	double **Resultado,**AuxResultado;
	
	/* ========================================
	 /*  Definicion de las variables de trabajo
	 ======================================== */
	
	// Ints
	I=*II;
	T=*TT;
	D=*DD;
	IT=I*T;
	
	// Matrices con datos
	YmA=array2srce(YA,I,T*D); //----------------- Liberar memoria
	YmB=array2srce(YB,I,T*D); //----------------- Liberar memoria
	SSi=array2srce(Si,IT,T); //------------------ Liberar memoria
	TamGrupos=array2srce(TamGrp, I, 1); //------- Liberar memoria
	Alfa=array2srce(alfaI, I, T); //------------- Liberar memoria
	
	// Matrices incializadas en 0's
	Mu=array2(T,1); //--------------------------- Liberar memoria
	AuxMu=array2(T,1); //------------------------ Liberar memoria
	AuxMu2=array2(T,1); //----------------------- Liberar memoria
	
	AuxMuFinal=array2(T,T); //------------------- Liberar memoria
	AuxMuFinal2=array2(T,1); //------------------ Liberar memoria
	AuxMuFinal3=array1(T); //------------------ Liberar memoria
	
	ind=arrayInt1(T); //------------------------- Liberar memoria
	MuA=array2(T,1);  //------------------------- Liberar memoria
	MuB=array2(T,1);  //------------------------- Liberar memoria
	
	SubYmA=array2(I,T);  //---------------------- Liberar memoria
	SubYmB=array2(I,T);  //---------------------- Liberar memoria
	SubYmAi=array2(T,1); //---------------------- Liberar memoria
	SubYmBi=array2(T,1); //---------------------- Liberar memoria
	rSubYmAi=array2(T,1); //--------------------- Liberar memoria
	rSubYmBi=array2(T,1); //--------------------- Liberar memoria
	
	Si1=array2(T,T); //-------------------------- Liberar memoria
	InvSi=array2(T,T); //------------------------ Liberar memoria
	Cest=array2(T,T); //------------------------- Liberar memoria
	InvCest=array2(T,T); //---------------------- Liberar memoria
	
	Resultado=array2(I,T); //------------------ Liberar memoria
	AuxResultado=array2(1,T); //----------------- Liberar memoria
	
	for(i=0;i<I;i++)
	{
		ki=TamGrupos[i][0];
		for(k=0;k<T;k++)
		{
			Mu[k][0]=0;
		}
		
		
		// Obtencion de los vectores MU
		for(j=0;j<D;j++)
		{
			// Definicion del vector de indices y de las medias de los grupos
			for(k=0;k<T;k++)
			{
				ind[k]=indices[k]+j;
				MuA[k][0]=0;
				MuB[k][0]=Alfa[i][k];
			}
			
			SubMat(YmA, I, ind, T,SubYmA);
			SubMat(YmB, I, ind, T,SubYmB);
			for(k=0;k<T;k++)
			{
				SubYmAi[k][0]=SubYmA[i][k];
				SubYmBi[k][0]=SubYmB[i][k];
			}
			MatRes(SubYmAi,MuA,T,1,rSubYmAi);
			MatRes(SubYmBi,MuB,T,1,rSubYmBi);
			
			MatSum(rSubYmAi, rSubYmBi, T, 1, AuxMu2);
			MatSum(Mu, AuxMu2, T, 1, AuxMu);
			
			for(k=0;k<T;k++)
			{
				Mu[k][0]=AuxMu[k][0];
				AuxMu[k][0]=0;
				AuxMu2[k][0]=0;
			}
		}
		
		a=T*i;
		//Obtencion de S.i.
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				Si1[j][k]=SSi[a+j][k];
			}
		}
		MatInv(Si1, TT, InvSi);
		
		// Multiplicacion de la matriz Cest
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				Cest[j][k]=2*D*ki*InvSi[j][k];
			}
		}
		MatInv(Cest, TT, InvCest);
		
		MatMult(InvCest,InvSi,T,T,T,AuxMuFinal);
		MatMult(AuxMuFinal, Mu, T, T, 1, AuxMuFinal2);
		
		for(j=0;j<T;j++)
		{
			AuxMuFinal3[j]=ki*AuxMuFinal2[j][0];
		}
		
		Crmnorm(1, T, AuxMuFinal3, InvCest, AuxResultado);
		
		for(k=0;k<T;k++)
			Resultado[i][k]=AuxResultado[0][k];
		
		//printf("%lf \n",ki);
	}
	
	contador=0;
	for(k=0;k<T;k++)
	{
		for(j=0;j<I;j++)
		{
			Resul[contador]=Resultado[j][k];
			contador=contador+1;
		}
	}
	
}

//===============================================================
//========    NP Alfa 0    =========================================
//===============================================================

void conditionalNPAlfa0(int *II,int *TT,int *DD,int *indices,double *TamGrp,double *YB,double *betaI,double *alfaI,double *Si,double *Resul)
{
	int I,D,T,i,j,k,a,contador;
	int *ind;
	double ki,varianza,M,media;
	double *unos;
	double **YmB,**SSi,**TamGrupos,**Beta,**Alfa;
	double **Si1,**InvSi,**UnosVec,**TUnosVec,**auxVar,**auxVar2,**MuB,**SubYmB;
	double **SubYmBi,**rSubYmBi,**mSubYmBi,**MB,**Resultado;
	
	/* ========================================
	 /*  Definicion de las variables de trabajo
	 ======================================== */
	
	// Ints
	I=*II;
	T=*TT;
	D=*DD;
	
	// Matrices con datos
	YmB=array2srce(YB,I,T*D); //----------------- Liberar memoria
	SSi=array2srce(Si,I*T,T); //------------------ Liberar memoria
	TamGrupos=array2srce(TamGrp, I, 1); //------- Liberar memoria
	Beta=array2srce(betaI, I, T); //------------- Liberar memoria
	Alfa=array2srce(alfaI, I, T); //----------- Liberar memoria
	
	// Matrices incializadas en 0's
	Si1=array2(T,T); //-------------------------- Liberar memoria
	InvSi=array2(T,T); //-------------------------- Liberar memoria
	UnosVec=array2(T,1); //-------------------------- Liberar memoria
	TUnosVec=array2(1,T); //-------------------------- Liberar memoria
	auxVar=array2(T,1); //-------------------------- Liberar memoria
	auxVar2=array2(1,1); //-------------------------- Liberar memoria
	
	ind=arrayInt1(T); //------------------------- Liberar memoria
	unos=array1(T); //-------------------------- Liberar memoria
	MuB=array2(T,1);  //------------------------- Liberar memoria
	
	SubYmB=array2(I,T);  //---------------------- Liberar memoria
	SubYmBi=array2(T,1); //---------------------- Liberar memoria
	rSubYmBi=array2(T,1); //--------------------- Liberar memoria
	mSubYmBi=array2(T,1); //--------------------- Liberar memoria
	MB=array2(1,1); //--------------------- Liberar memoria
	
	Resultado=array2(I,1); //------------------ Liberar memoria
	
	for(j=0;j<T;j++)
	{
		unos[j]=1;
		UnosVec[j][0]=1;
		TUnosVec[0][j]=1;
	}
	
	for(i=0;i<I;i++)
	{
		ki=TamGrupos[i][0];
		a=T*i;
		
		// Definicion de la matriz (1/ki)S.i y su inversa
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
				Si1[j][k]=(1/ki)*SSi[a+j][k];
		}
		MatInv(Si1, TT, InvSi);
		MatMult(InvSi, UnosVec, T, T, 1, auxVar);
		MatMult(TUnosVec, auxVar, 1, T, 1, auxVar2);
		
		M=0;
		varianza=0;
		// Ciclo dentro de dosis
		for(j=0;j<D;j++)
		{
			varianza=varianza+auxVar2[0][0];
			
			// Definicion del vector de indices y de las medias de los grupos
			for(k=0;k<T;k++)
			{
				ind[k]=indices[k]+j;
				MuB[k][0]=Beta[i][k]+Alfa[i][k];
			}
			SubMat(YmB, I, ind, T,SubYmB);
			for(k=0;k<T;k++)
			{
				SubYmBi[k][0]=SubYmB[i][k];
			}
			MatRes(SubYmBi,MuB,T,1,rSubYmBi);
			
			MatMult(InvSi, rSubYmBi, T, T, 1, mSubYmBi);
			
			MatMult(TUnosVec, mSubYmBi, 1, T, 1, MB);
			
			M=M+MB[0][0];
		}
		
		varianza=1/varianza;
		media=M*varianza;
		varianza=sqrt(varianza);
		
		GetRNGstate();	
		Resultado[i][0]=rnorm(media,varianza);
		PutRNGstate();
		//printf("%lf \n",Resultado[i][0]);
	}
	
	contador=0;
	for(k=0;k<I;k++)
	{
		Resul[contador]=Resultado[k][0];
		contador=contador+1;
	}
		
}

//===============================================================
//=========    NP Alfa    =========================================
//===============================================================

void conditionalNPAlfa(int *II,int *TT,int *DD,int *indices,double *TamGrp,double *YB,double *betaI,double *alfa0I,double *Si,double *Sai,double *rho,double *mm,double *Resul)
{
	int i,j,k,I,T,D,a,contador;
	int *ind;
	double ki,suma,auxSuma,zuma,auxZuma,MM,rzn_log,numI,demI,fracI,rzn,prob,alea;
	double *ceros,*auxrSubYmBi,*AuxMuFinal3;
	double **YmB,**SSi,**SSa,**TamGrupos,**Beta,**Alfa0,**Ro;
	double **SS_Si,**SS_Si_Bis,**SS_No,**SS_No_Bis,**Mu,**AuxMu,**MuB,**Si1,**Si2,**InvSi,**InvSa,**Cest,**InvCest;
	double **SubYmB,**SubYmBi,**rSubYmBi,**AuxMuFinal,**AuxMuFinal2,**Resultado,**AuxResultado;
	
	/* ========================================
	 /*  Definicion de las variables de trabajo
	 ======================================== */
	
	// Ints
	I=*II;
	T=*TT;
	D=*DD;
	MM=*mm;
	
	// Matrices con datos
	YmB=array2srce(YB,I,T*D); //----------------- Liberar memoria
	SSi=array2srce(Si, I*T, T); //--------------- Liberar memoria
	SSa=array2srce(Sai,T,T); //------------------ Liberar memoria
	TamGrupos=array2srce(TamGrp, I, 1); //------- Liberar memoria
	Beta=array2srce(betaI, I, T); //------------- Liberar memoria
	Alfa0=array2srce(alfa0I, I, T); //------------- Liberar memoria
	Ro=array2srce(rho, I, 1); //----------------- Liberar memoria
	
	
	// Matrices incializadas en 0's
	SS_Si=array2(T,T); //------------------------ Liberar memoria
	SS_Si_Bis=array2(T,T); //-------------------- Liberar memoria
	SS_No=array2(T,T); //------------------------ Liberar memoria
	SS_No_Bis=array2(T,T); //-------------------- Liberar memoria
	Mu=array2(T,1); //--------------------------- Liberar memoria
	AuxMu=array2(T,1); //------------------------ Liberar memoria 
	MuB=array2(T,1);  //------------------------- Liberar memoria
	Si1=array2(T,T); //-------------------------- Liberar memoria
	Si2=array2(T,T); //-------------------------- Liberar memoria
	InvSi=array2(T,T); //------------------------ Liberar memoria
	InvSa=array2(T,T); //------------------------ Liberar memoria
	Cest=array2(T,T); //------------------------- Liberar memoria
	InvCest=array2(T,T); //---------------------- Liberar memoria
	
	SubYmB=array2(I,T); //----------------------- Liberar memoria
	SubYmBi=array2(T,1); //---------------------- Liberar memoria
	rSubYmBi=array2(T,1); //--------------------- Liberar memoria
	auxrSubYmBi=array1(T); //-------------------- Liberar memoria
	
	AuxMuFinal=array2(T,T); //------------------- Liberar memoria
	AuxMuFinal2=array2(T,1); //------------------ Liberar memoria
	AuxMuFinal3=array1(T); //-------------------- Liberar memoria
	
	Resultado=array2(I,T); //-------------------- Liberar memoria
	AuxResultado=array2(1,T); //----------------- Liberar memoria
	
	
	ind=arrayInt1(T); //------------------------- Liberar memoria
	ceros=array1(T); //-------------------------- Liberar memoria
	
	for(i=0;i<T;i++)
	{
		ceros[i]=0;
	}
	
	for(i=0;i<I;i++)
	{
		a=T*i;
		ki=TamGrupos[i][0];
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				SS_Si[j][k]=(SSi[a+j][k]/ki)+SSa[j][k];
				SS_Si_Bis[j][k]=(SSi[a+j][k]/ki)+SSa[j][k];
				SS_No[j][k]=(SSi[a+j][k]/ki);
				SS_No_Bis[j][k]=(SSi[a+j][k]/ki);
				Si1[j][k]=SSi[a+j][k];
				Si2[j][k]=SSa[j][k];
			}
			Mu[j][0]=0;
		}
		
		suma=0;
		zuma=0;
		// Obtencion de los vectores MU
		for(j=0;j<D;j++)
		{
			// Definicion del vector de indices y de las medias de los grupos
			for(k=0;k<T;k++)
			{
				ind[k]=indices[k]+j;
				MuB[k][0]=Beta[i][k]+Alfa0[i][k];
			}
			
			SubMat(YmB, I, ind, T,SubYmB);
			for(k=0;k<T;k++)
				SubYmBi[k][0]=SubYmB[i][k];
			MatRes(SubYmBi,MuB,T,1,rSubYmBi);
			
			for(k=0;k<T;k++)
				auxrSubYmBi[k]=rSubYmBi[k][0];
			
			MatSum(Mu, rSubYmBi, T, 1, AuxMu);
			
			auxSuma=Cdmnorm(TT, auxrSubYmBi, ceros, SS_Si, SS_Si_Bis);
			auxZuma=Cdmnorm(TT, auxrSubYmBi, ceros, SS_No, SS_No_Bis);
			
			for(k=0;k<T;k++)
			{
				Mu[k][0]=AuxMu[k][0];
				AuxMu[k][0]=0;
			}
			
			suma=suma+auxSuma;
			zuma=zuma+auxZuma;
		}
		
		numI=Ro[i][0]*MM;
		demI=1-(Ro[i][0]*MM);
		fracI=numI/demI;
		rzn_log=suma-zuma+log(fracI);
		
		if(rzn_log>500)
		{
			prob=1;
		}else if(rzn_log<(-500))
		{
			prob=0;
		}else{
			rzn=exp(rzn_log);
			prob=rzn/(1+rzn);
		}
		
		MatInv(Si1, TT, InvSi);
		MatInv(Si2, TT, InvSa);
		// Multiplicacion de la matriz Cest
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				Cest[j][k]=D*ki*InvSi[j][k]+InvSa[j][k];
			}
		}
		MatInv(Cest, TT, InvCest);
		
		MatMult(InvCest,InvSi,T,T,T,AuxMuFinal);
		MatMult(AuxMuFinal, Mu, T, T, 1, AuxMuFinal2);
		
		for(j=0;j<T;j++)
		{
			AuxMuFinal3[j]=ki*AuxMuFinal2[j][0];
		}
		
		GetRNGstate();
		alea=runif(0,1);
		PutRNGstate();
		
		if(alea<prob)
		{
			Crmnorm(1, T, AuxMuFinal3, InvCest, AuxResultado);
			for(k=0;k<T;k++)
				Resultado[i][k]=AuxResultado[0][k];
		}else{
			for(k=0;k<T;k++)
				Resultado[i][k]=0;
		}
	}
	
	contador=0;
	for(k=0;k<T;k++)
	{
		for(j=0;j<I;j++)
		{
			Resul[contador]=Resultado[j][k];
			contador=contador+1;
		}
	}
	
}

//===============================================================
//=========    Alfa    =========================================
//===============================================================

void conditionalALLAlfa(int *II,int *TT,int *DD,int *indices,double *TamGrp,double *YB,double *betaI,double *Si,double *Sai,double *rho,double *mm,double *Resul)
{
	int i,j,k,I,T,D,a,contador;
	int *ind;
	double ki,suma,auxSuma,zuma,auxZuma,MM,rzn_log,numI,demI,fracI,rzn,prob,alea;
	double *ceros,*auxrSubYmBi,*AuxMuFinal3;
	double **YmB,**SSi,**SSa,**TamGrupos,**Beta,**Ro;
	double **SS_Si,**SS_Si_Bis,**SS_No,**SS_No_Bis,**Mu,**AuxMu,**MuB,**Si1,**Si2,**InvSi,**InvSa,**Cest,**InvCest;
	double **SubYmB,**SubYmBi,**rSubYmBi,**AuxMuFinal,**AuxMuFinal2,**Resultado,**AuxResultado;
	
	/* ========================================
	 /*  Definicion de las variables de trabajo
	 ======================================== */
	
	// Ints
	I=*II;
	T=*TT;
	D=*DD;
	MM=*mm;
	
	// Matrices con datos
	YmB=array2srce(YB,I,T*D); //----------------- Liberar memoria
	SSi=array2srce(Si, I*T, T); //--------------- Liberar memoria
	SSa=array2srce(Sai,T,T); //------------------ Liberar memoria
	TamGrupos=array2srce(TamGrp, I, 1); //------- Liberar memoria
	Beta=array2srce(betaI, I, T); //------------- Liberar memoria
	Ro=array2srce(rho, I, 1); //----------------- Liberar memoria
	
	
	// Matrices incializadas en 0's
	SS_Si=array2(T,T); //------------------------ Liberar memoria
	SS_Si_Bis=array2(T,T); //-------------------- Liberar memoria
	SS_No=array2(T,T); //------------------------ Liberar memoria
	SS_No_Bis=array2(T,T); //-------------------- Liberar memoria
	Mu=array2(T,1); //--------------------------- Liberar memoria
	AuxMu=array2(T,1); //------------------------ Liberar memoria 
	MuB=array2(T,1);  //------------------------- Liberar memoria
	Si1=array2(T,T); //-------------------------- Liberar memoria
	Si2=array2(T,T); //-------------------------- Liberar memoria
	InvSi=array2(T,T); //------------------------ Liberar memoria
	InvSa=array2(T,T); //------------------------ Liberar memoria
	Cest=array2(T,T); //------------------------- Liberar memoria
	InvCest=array2(T,T); //---------------------- Liberar memoria
	
	SubYmB=array2(I,T); //----------------------- Liberar memoria
	SubYmBi=array2(T,1); //---------------------- Liberar memoria
	rSubYmBi=array2(T,1); //--------------------- Liberar memoria
	auxrSubYmBi=array1(T); //-------------------- Liberar memoria
	
	AuxMuFinal=array2(T,T); //------------------- Liberar memoria
	AuxMuFinal2=array2(T,1); //------------------ Liberar memoria
	AuxMuFinal3=array1(T); //-------------------- Liberar memoria
	
	Resultado=array2(I,T); //-------------------- Liberar memoria
	AuxResultado=array2(1,T); //----------------- Liberar memoria
	
	
	ind=arrayInt1(T); //------------------------- Liberar memoria
	ceros=array1(T); //-------------------------- Liberar memoria
	
	for(i=0;i<T;i++)
	{
		ceros[i]=0;
	}
	
	for(i=0;i<I;i++)
	{
		a=T*i;
		ki=TamGrupos[i][0];
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				SS_Si[j][k]=(SSi[a+j][k]/ki)+SSa[j][k];
				SS_Si_Bis[j][k]=(SSi[a+j][k]/ki)+SSa[j][k];
				SS_No[j][k]=(SSi[a+j][k]/ki);
				SS_No_Bis[j][k]=(SSi[a+j][k]/ki);
				Si1[j][k]=SSi[a+j][k];
				Si2[j][k]=SSa[j][k];
			}
			Mu[j][0]=0;
		}
		
		suma=0;
		zuma=0;
		// Obtencion de los vectores MU
		for(j=0;j<D;j++)
		{
			// Definicion del vector de indices y de las medias de los grupos
			for(k=0;k<T;k++)
			{
				ind[k]=indices[k]+j;
				MuB[k][0]=Beta[i][k];
			}
			
			SubMat(YmB, I, ind, T,SubYmB);
			for(k=0;k<T;k++)
				SubYmBi[k][0]=SubYmB[i][k];
			MatRes(SubYmBi,MuB,T,1,rSubYmBi);
			
			for(k=0;k<T;k++)
				auxrSubYmBi[k]=rSubYmBi[k][0];
			
			MatSum(Mu, rSubYmBi, T, 1, AuxMu);
			
			auxSuma=Cdmnorm(TT, auxrSubYmBi, ceros, SS_Si, SS_Si_Bis);
			auxZuma=Cdmnorm(TT, auxrSubYmBi, ceros, SS_No, SS_No_Bis);
			
			for(k=0;k<T;k++)
			{
				Mu[k][0]=AuxMu[k][0];
				AuxMu[k][0]=0;
			}
			
			suma=suma+auxSuma;
			zuma=zuma+auxZuma;
		}
		
		numI=Ro[i][0]*MM;
		demI=1-(Ro[i][0]*MM);
		fracI=numI/demI;
		rzn_log=suma-zuma+log(fracI);
		
		if(rzn_log>500)
		{
			prob=1;
		}else if(rzn_log<(-500))
		{
			prob=0;
		}else{
			rzn=exp(rzn_log);
			prob=rzn/(1+rzn);
		}
		
		MatInv(Si1, TT, InvSi);
		MatInv(Si2, TT, InvSa);
		// Multiplicacion de la matriz Cest
		for(j=0;j<T;j++)
		{
			for(k=0;k<T;k++)
			{
				Cest[j][k]=D*ki*InvSi[j][k]+InvSa[j][k];
			}
		}
		MatInv(Cest, TT, InvCest);
		
		MatMult(InvCest,InvSi,T,T,T,AuxMuFinal);
		MatMult(AuxMuFinal, Mu, T, T, 1, AuxMuFinal2);
		
		for(j=0;j<T;j++)
		{
			AuxMuFinal3[j]=ki*AuxMuFinal2[j][0];
		}
		
		GetRNGstate();
		alea=runif(0,1);
		PutRNGstate();
		
		if(alea<prob)
		{
			Crmnorm(1, T, AuxMuFinal3, InvCest, AuxResultado);
			for(k=0;k<T;k++)
				Resultado[i][k]=AuxResultado[0][k];
		}else{
			for(k=0;k<T;k++)
				Resultado[i][k]=0;
		}
	}
	
	contador=0;
	for(k=0;k<T;k++)
	{
		for(j=0;j<I;j++)
		{
			Resul[contador]=Resultado[j][k];
			contador=contador+1;
		}
	}
}
