#include    <math.h>
#include    <stdio.h>
#include    <stdlib.h>

double	Likelihood(int	T[][4]){
	int i,	j;
	int c[2]={0,	0};
	double  l=0;


	for(i=0;    i<4;    i++){
		for(j=0;    j<2;    j++){
	    	c[j]+=T[j][i];
			l-=lgamma(T[j][i]+1);
		}
        l+=lgamma(T[0][i]+T[1][i]+1);
	}
	l+=lgamma(c[0]+1)+lgamma(c[1]+1)-lgamma(c[0]+c[1]+1);
	return  exp(l);
}

double	ExactTest42(int	T[][4]){
	double	l,	m,	p=0;
	int	x,	y,	z,	i,	j;
	int X[2]={0,	0},	Y[4]={0,	0,	0,	0};
	int t[2][4];

	for(i=0;    i<4;    i++)	for(j=0;    j<2;    j++){
		X[j]+=T[j][i];
		Y[i]+=T[j][i];
	}

	m=Likelihood(T);
	for(x=0;	x<=Y[0];	x++)
	for(y=0;    y<=Y[1]&&	X[0]-x-y>=0;	y++)
    for(z=0;    z<=Y[2]&&	X[0]-x-y-z>=0;	z++){
		t[0][0]=x;  t[0][1]=y;  t[0][2]=z;  t[0][3]=X[0]-x-y-z;
		t[1][0]=Y[0]-t[0][0];
		t[1][1]=Y[1]-t[0][1];
		t[1][2]=Y[2]-t[0][2];
		t[1][3]=Y[3]-t[0][3];

  		l=Likelihood(t);
  		if(l<=m)	p+=l;
	}
	return	p;
}
/*
int main(void){
	int T[2][4]={	29,	88,	117,	734,
					0,	99,	131,	608};
	printf("%g	%g\n",  Likelihood(T),	ExactTest42(T));
	system("pause");
}
*/
