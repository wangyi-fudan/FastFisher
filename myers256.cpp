#include	<x86intrin.h>
#include	<iostream>
#include	<stdint.h>
using	namespace	std;
typedef	long	long	int	v4di	__attribute__	((__vector_size__	(32)));
union	YMM {	uint64_t	x[4];	__uint128_t	y[2];	v4di	z;	};

   
struct	Myers256{
	YMM	mask[5];
	uint16_t	len,	wor,	shift;
	inline	void	ShiftLeft(YMM	&X){
	   X.z=_mm256_or_si256(_mm256_slli_epi64(X.z,	1),	_mm256_blend_epi32(_mm256_setzero_si256(),	_mm256_permute4x64_epi64(_mm256_srli_epi64(X.z, 63), 0x93),	0xFC));
	}
	void	encode(const	void	*S,	uint64_t	L){
		len=L;	wor=(len-1)>>6ULL;	shift=(len-1)&63;
		for(uint64_t	i=0ULL;	i<5ULL;	i++)	mask[i].x[0]=mask[i].x[1]=mask[i].x[2]=mask[i].x[3]=0ULL;
		for(uint64_t	i=0ULL;	i<len;	i++)	mask[((const	uint8_t	*)S)[i]].x[i>>6ULL]|=1ULL<<(i&63ULL);
	}
	uint16_t	distance(const	unsigned	*S,	unsigned	B,	unsigned	E,	unsigned	&P){
		YMM	pv={~0ULL,~0ULL,~0ULL,~0ULL},	ph,	mv={0ULL,0ULL,0ULL,0ULL},	mh,	xv,	xh,	*eq,	t1=pv;
		uint16_t	s=len,	b=len+1;
		for(unsigned	i=B;	i<E;	i++){
			eq=mask+((S[i>>4]>>((i&15U)<<1U))&3U);
			xv.z=__builtin_ia32_por256(eq->z,	mv.z);
			xh.z=__builtin_ia32_andsi256(eq->z,	pv.z);
			xh.y[0]+=pv.y[0];	xh.y[1]+=pv.y[1]+(xh.y[0]<pv.y[0]);
			xh.z=__builtin_ia32_pxor256(xh.z,	pv.z);
			xh.z=__builtin_ia32_por256(xh.z,	eq->z);
			ph.z=__builtin_ia32_por256(xh.z,	pv.z);
			ph.z=__builtin_ia32_pxor256(t1.z,	ph.z);
			ph.z=__builtin_ia32_por256(mv.z,	ph.z);
			mh.z=__builtin_ia32_andsi256(pv.z,	xh.z);
			s+=(ph.x[wor]>>shift)&1ULL;	s-=(mh.x[wor]>>shift)&1ULL;
			ShiftLeft(ph);	ShiftLeft(mh);
			xv.z=__builtin_ia32_por256(eq->z,	mv.z);
			pv.z=__builtin_ia32_por256(xv.z,	ph.z);
			pv.z=__builtin_ia32_pxor256(t1.z,	pv.z);
			pv.z=__builtin_ia32_por256(mh.z,	pv.z);
			mv.z=__builtin_ia32_andsi256(ph.z,	xv.z);
			if(s<b){	b=s;	P=i;	}
		}
		return	b;
	}		
};

