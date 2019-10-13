#ifndef	MST_Included
#define MST_Included
#include    <cfloat>
double	MST(float	*D,	int	N){
	int	*V=new	int[N];	float	*M=new	float[N];
	int	i,	j,	k;	float    m;	double	s=0;
	k=0;	V[0]=-1;	for(j=1;	j<N;	j++){	V[j]=0;	M[j]=D[0*N+j];	}
	for(i=1;	i<N;	i++){
		m=FLT_MAX;	
		for(j=0;	j<N;	j++)	if(M[j]<=m&&V[j]>=0){	m=M[j];	k=j;	}
		s+=D[k*N+V[k]];	V[k]=-1;
		for (j=0;	j<N;	j++)	if(V[j]>=0&&D[k*N+j]<=M[j]){	V[j]=k;	M[k]=D[k*N+j];	}
	}
	delete	[]	V;	delete	[]	M;
	return	s;
}
#endif

