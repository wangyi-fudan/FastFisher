#include	<algorithm>
#include	<cmath>
using	namespace	std;

struct	Sort{
	double	x;
	size_t	i;	
	bool	operator()(Sort	X,	Sort	Y){	return	X.x<Y.x;	}
};

static	inline	void	Transform2CDF(const	double	*x0,	double	*x1,	size_t	n){
	Sort	*s=new	Sort[n];
	for(size_t	i=0;	i<n;	i++){	s[i].x=x0[i];	s[i].i=i;	}
	sort(s,s+n,Sort());
	for(size_t	i=0;	i<n;){
		size_t	j;
		for(j=i;	j<n&&s[i].x==s[j].x;	j++);
		double	r=0.5*(i+j)/n;
		for(size_t	k=i;	k<j;	k++)	x1[s[k].i]=r;
		i=j;
	}
	delete	[]	s;
}

static	inline	double	GeneralizedCorrelationLRT(const	double	*x,	const	double	*y,	size_t	n){
	double	*cx=new	double	[n],	*cy=new	double	[n];
	Transform2CDF(x,cx,n);	Transform2CDF(y,cy,n);
	double	tau=1.0/sqrt(n),	l=0;
	for(size_t	i=0;	i<n;	i++){
		size_t	sxy=0,	sx=0,	sy=0;
		for(size_t	j=0;	j<n;	j++){
			sx+=(fabs(cx[i]-cx[j])<tau);
			sy+=(fabs(cy[i]-cy[j])<tau);
			sxy+=(fabs(cx[i]-cx[j])<tau)&&(fabs(cy[i]-cy[j])<tau);
		}
		l+=log((double)sxy*n/sx/sy);
	}
	delete	[]	cx;	delete	[]	cy;
	return	2*l;
}

double	GeneralizedCorrelation(const	double	*x,	const	double	*y,	size_t	n){
	double	chi2=GeneralizedCorrelationLRT(x,y,n);
	double	*sx=new	double	[n],	*sy=new	double	[n];
	for(size_t	i=0;	i<n;	i++){	sx[i]=x[i];	sy[i]=y[i];	}
	sort(sx,sx+n);	sort(sy,sy+n);
	double	chi2max=GeneralizedCorrelationLRT(sx,sy,n);
	delete	[]	sx;	delete	[]	sy;
	return	chi2/chi2max;
}

#include	<iostream>
int	main(void){
	double	x[10]={1,1,1,2,2,2,2,3,3,3};
	double	y[10]={1,7,3,4,5,6,7,8,9,1};
	cerr<<GeneralizedCorrelationLRT(x,y,10)<<'\n';
}

