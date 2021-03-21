#include	<algorithm>
#include	<iostream>
#include	<stdlib.h>
#include	"wyhash.h"
#include	<cmath>
using	namespace	std;

struct	FastSpearmanSort{
	double	v;
	int64_t	i;
	bool	operator()(FastSpearmanSort	X,	FastSpearmanSort	Y){	return	X.v<Y.v;	}
};

static	inline	void	FastSpearmanRank(const	double	*x,	int64_t	*r,	unsigned	n){
	struct	FastSpearmanSort*	vs=(struct  FastSpearmanSort*)malloc(n*sizeof(struct  FastSpearmanSort));
	for(unsigned	i=0;	i<n;	i++){	vs[i].v=x[i];	vs[i].i=i;	}
	sort(vs, vs+n, FastSpearmanSort());
	for(unsigned	i=0;	i<n;){
		unsigned	j;
		for(j=i+1;	j<n&&vs[i].v==vs[j].v;	j++);
		unsigned	rank=i+j;
		for(unsigned	k=i;	k<j;	k++)	r[vs[k].i]=rank;
		i=j;
	}
	free(vs);
}

double	FastSpearman(const	double	*x,	const	double	*y,	unsigned	n,	uint64_t	N,	double	*rho){
	int64_t	*a=(int64_t*)malloc(n*sizeof(int64_t)),	*b=(int64_t*)malloc(n*sizeof(int64_t));
	int64_t	obs=0, per=0, ma=0,mi=0;;
	FastSpearmanRank(x,a,n);	FastSpearmanRank(y,b,n);
	double	sx=0,	sxx=0,	sy=0,	syy=0,	sxy=0;
	for(unsigned	i=0;	i<n;	i++){
		obs+=a[i]*b[i];	sx+=a[i];	sxx+=a[i]*a[i];	sy+=b[i];	syy+=b[i]*b[i];	sxy+=a[i]*b[i];
	}
	sx/=n;	sy/=n;	sxx=sxx/n-sx*sx;	syy=syy/n-sy*sy;
	*rho=sxx>0&&syy>0?(sxy/n-sx*sy)/sqrt(sxx*syy):0;

	sort(a,a+n);	sort(b,b+n);
	for(unsigned	i=0;	i<n;	i++)	ma+=a[i]*b[i];
	sort(a,a+n,greater<int64_t>());	
	for(unsigned	i=0;	i<n;	i++)	mi+=a[i]*b[i];
	double	beta=logf(0.1*N)/max(ma-obs,obs-mi);

	uint64_t	seed=0;
	for(unsigned	i=n-1;	i;	i--)	swap(a[i],a[wy2u0k(wyrand(&seed),i+1)]);
	for(unsigned	i=0;	i<n;	i++)	per+=a[i]*b[i];

	double	ge[2]={},	le[2]={};
	for(uint64_t	p=0;	p<N;	p++){
		unsigned	i,	j;
		do{	uint64_t	r=wyrand(&seed);	i=((r>>32)*n)>>32;	j=((r&0xffffffff)*n)>>32;	}while(i==j);
		int64_t	delta=(a[j]-a[i])*(b[i]-b[j]),	d0=abs(per-obs),	d1=abs(per+delta-obs);
		if(wy2u01(wyrand(&seed))<expf(beta*(d0-d1))){	per+=delta;	swap(a[i],	a[j]);	}
		if(p<N/10)	continue;
		double	prob=expf(beta*abs(per-obs));
		ge[per>=obs]+=prob;	le[per<=obs]+=prob;
	}
	free(a);	free(b);
	return	2*min(ge[1]/(ge[1]+ge[0]),le[1]/(le[1]+le[0]));
}

int	main(void){
	double	x[11]={1,2,3,4,5,6,7,8,9,10,11},y[11]={1,2,3,9,5,6,7,8,4,10,11},	rho;
	cerr<<FastSpearman(x,y,11,100000ull,&rho)<<'\n';
	cerr<<rho<<'\n';
	return	0;
}
