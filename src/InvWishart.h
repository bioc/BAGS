/*
 *  InvWishart.h
 *
 *  Created by Alejandro Quiroz on 3/20/11.
 *
 */

#include "Densidades.h"


void Criwish(int v, int *pp, double **Sg, double **ICrw)
{
	int p;
	double **rSg, **Crw;
	
	p=*pp;
	rSg=array2(p,p);
	Crw=array2(p,p);
	
	MatInv(Sg,pp,rSg);	
	Crwish(v,p,rSg,Crw);
	MatInv(Crw,pp,ICrw);
	
}
