/*
  Name:	HWE(exact test)
  Copyright:	Forever
  Author:	Fu wenqing
  Date: 10-03-06
  Description:	This function is used to calculate HWE Exact Test
*/

#ifndef	ExactTestHWE_Included
#define	ExactTestHWE_Included

#include <math.h>

double Pr(int a,int b,int c){
       int n,x,y,i;
       double pi;
       n=a+b+c;
       x=2*a+b;
       y=2*c+b;
       pi=lgamma(n+1)+lgamma(x+1)+lgamma(y+1)-lgamma(a+1)-lgamma(b+1)-lgamma(c+1)-lgamma(2*n+1);
	   for(i=1;i<=b;i++) pi=pi+log(2);
	   pi=exp(pi);
       return(pi);
}

double HWE(int a,int b,int c){
       int r,s,t,u,v;
       double p,fo,f;
       r=2*a+b;
       s=2*c+b;
       t=a;u=b;v=c;
       fo=Pr(a,b,c);
       p=fo;
       b=b-2;a=a+1;c=c+1;
       while(b>=0){
           f=Pr(a,b,c);
           b=b-2;a=a+1;c=c+1;
           if(f<=fo) p+=f;
       }
       u=u+2;t=t-1;v=v-1;
       while((t>=0)&&(v>=0)){
           f=Pr(t,u,v);
           u=u+2;t=t-1;v=v-1;
           if(f<=fo) p+=f;
       }
       return(p);
}

#include	<iostream>
using	namespace	std;
int	main(void){
	int	AA,	Aa,	aa;
	cout<<"Welcome to Two Allele HWE Exact Test!\n";
	cout<<"Wenqing Fu and Yi Wang @ Fudan University.\n";
	cout<<"Enter count of AA, Aa, aa:\t";
	cin>>AA>>Aa>>aa;
	cout<<HWE(AA,	Aa,	aa)<<endl;
//	system("pause");	
}

#endif
