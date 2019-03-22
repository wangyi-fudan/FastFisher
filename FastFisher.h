/*Author:	Wang Yi*/
#ifndef	FastFisher_Included
#define	FastFisher_Included
#include	<math.h>
long	double	FastFisher(unsigned	AA,	unsigned	AB,	unsigned	BA,	unsigned	BB){
	long	double	l0=expl(lgammal(AA+AB+1)+lgammal(BA+BB+1)+lgammal(AA+BA+1)+lgammal(AB+BB+1)-lgammal(AA+AB+BA+BB+1)-lgammal(AA+1)-lgammal(AB+1)-lgammal(BA+1)-lgammal(BB+1));
	long	double	pval=l0,	cut=l0*1.000000000001L,	l,	aa,	ab,	ba,	bb;
	long	double	beg=AA>BB?AA-BB:0,	end=AB<BA?AA+AB:AA+BA;
	for(aa=AA,	ab=AB,	ba=BA,	bb=BB,	l=l0;	aa>beg;	aa-=1,	ab+=1,	ba+=1,	bb-=1){	
		l*=aa*bb/((ab+1)*(ba+1));	if(l<cut)	pval+=l;
	}
	for(aa=AA,	ab=AB,	ba=BA,	bb=BB,	l=l0;	aa<end;	aa+=1,	ab-=1,	ba-=1,	bb+=1){	
		l*=ab*ba/((aa+1)*(bb+1));	if(l<cut)	pval+=l;	
	}
	return	pval;
}
#endif
