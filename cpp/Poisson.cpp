#include "Poisson.h"

float PI = 3.1415927;
float LZERO = -1e10;
float LSMALL = LZERO/2;
float SMALL = exp(LSMALL);
float minLogExp = -log(-LZERO);

float LogPoissonTailProb(float n, float lambda){
	float logprob = LZERO;
	float plogprob;
	do{
		plogprob = logprob;
		logprob = LAdd(plogprob, LogPoissonPDF(n++,lambda));
	} while(logprob-plogprob < 0.01);
	return logprob;
}

float LogPoissonPDF(float k, float lambda){
	float logk_factorial = k==0 ? 0 : k * (log (k) ) - k + 0.5 * (log (2*PI*k) ) ;
	float log_lambda = lambda<=0 ? LSMALL:log(lambda);
	float logp = k*log_lambda - lambda - logk_factorial;
	return logp;
}

float LAdd(float x, float y){
	float temp, diff, z;
	if(x < y){
		temp = x;
		x = y;
		y = temp;
	}
	diff = y - x;
	if(diff < minLogExp)
		return x<LSMALL?LZERO:x;
	else{
		z = exp(diff);
		return x+log(1.0+z);
	}
	return z;
}
