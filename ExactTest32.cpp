/*
  Name:	Exact Test
  Copyright:	Forever
  Author:	Wangyi
  Date: 18-04-05 22:01
  Description:	This function is used to calculate Fisher' Exact Test
*/

#ifndef	ExactTest32_Included
#define	ExactTest32_Included

#include	<math.h>

double	Likelihood32(int	A0,	int	B0,	int	C0,	int	A1,	int	B1,	int	C1,	double	M){
	return	exp(M-lgamma(A0+1)-lgamma(B0+1)-lgamma(C0+1)-lgamma(A1+1)-lgamma(B1+1)-lgamma(C1+1));
}

double	ExactTest32(int	A0,	int	B0,	int	C0,	int	A1,	int	B1,	int	C1){
	double	l,	m,	M,	p;
	int	a0,	b0,	c0,	a1,	b1,	c1;

	M=lgamma(A0+B0+C0+1)+lgamma(A1+B1+C1+1)+lgamma(A0+A1+1)+lgamma(B0+B1+1)+lgamma(C0+C1+1)-lgamma(A0+B0+C0+A1+B1+C1+1);

	p=m=Likelihood32(A0,	B0,	C0,	A1,	B1,	C1,	M);
	for(a0=0;   a0<=A0+A1;  a0++)   for(b0=0;   b0<=B0+B1;   b0++)	if(a0!=A0||b0!=B0){
		c0=A0+B0+C0-a0-b0;
		a1=A0+A1-a0;
		b1=B0+B1-b0;
		c1=A1+B1+C1-a1-b1;
		if(c0>=0&&c1>=0){
  			l=Likelihood32(a0,	b0,	c0,	a1,	b1,	c1,	M);
  			if(l<=m)	p+=l;
		}
	}
	return	p;
}

#include	<stdio.h>
#include	<stdlib.h>
int	main(void){
	int a0,	b0,	c0,	a1,	b1,	c1;
	printf("3X2列表Fisher Exact Test:\na0	b0	c0\na1	b1	c1\n");
	printf("请按照a0 b0 c0 a1 b1 c1依次输入你的3*2列表的六个值(用空格分开）：\n");
	scanf("%d%d%d%d%d%d",	&a0,	&b0,	&c0,	&a1,	&b1,	&c1);
	printf("你的表格是:\n%d	%d	%d\n%d	%d	%d\n\n",	a0,	b0,	c0,	a1,	b1,	c1);
	printf("P_Value=%g\n",	ExactTest32(a0,	b0,	c0,	a1,	b1,	c1));
	system("pause");
}

#endif

