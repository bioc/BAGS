//
//  libmat.h
//  
//
//  Created by Alejandro Quiroz on 12/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

/*
 *  libmat.h
 *
 *  Created by Alejandro Quiroz on 3/20/11.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

int *arrayInt1(int row);
double *array1(int row);
double **array2(int row, int col);
double **array2srce(double *srce, int row, int col);
double prodArray1(double *v1, double *v2, int m);
void array2toArray1(double **a, double *v, int nrow, int cols);
double *array2toVec(double **a, double *v ,int nrow, int cols);

/***************************************************************/

int *arrayInt1(int row)
{
	int *mat;
	
	mat=(int*)R_alloc(row,sizeof(int));
	return mat;
}


double *array1(int row)
{
	double *mat;
	
	mat=(double*)R_alloc(row,sizeof(double));
	return mat;
}

double **array2(int row, int col)
{
	double *ptr1, **mat;
	int i;
	
	ptr1=(double*)R_alloc(row*col,sizeof(double));
	mat=(double**)R_alloc(row,sizeof(double*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	return mat;
}

double **array2srce(double *srce, int row, int col)
{
	double *ptr1, **mat;
	int i, j, k;
	
	ptr1=(double*)R_alloc(row*col,sizeof(double));
	mat=(double**)R_alloc(row,sizeof(double*));
	for (i=0; i<row; i++)
		mat[i]=ptr1+(i*col);
	
	k=0;
	for(j=0; j<col; j++)
		for(i=0; i<row; i++){
			mat[i][j]=srce[k];
			k++;
		}
	return mat;
}


double prodArray1(double *v1, double *v2, int m)
{
	double resul;
	int i;
	
	resul=0;
	for (i=0; i<m; i++)
		resul+=(v2[i]*v1[i]);
	return resul;
}

void array2toArray1(double **a, double *v, int nrow, int cols)
{
	int j;
	
	for(j=0; j<cols; j++)
		v[j]=a[nrow][j];
}

double *array2toVec(double **a, double *v ,int nrow, int cols)
{
	int i,j;
	int contador;
	
	contador=0;
	for(i=0;i<nrow;i++)
	{
    	for(j=0; j<cols; j++)
		{
			v[contador]=a[i][j];
			contador=contador+1;
		}
	}
	return v;
}
