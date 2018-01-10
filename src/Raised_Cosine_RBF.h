/*
 * Raised_Cosine_RBF.h
 *
 *  Created on: Jun 22, 2016
 *      Author: chaturvedi
 */

#ifndef SRC_RAISED_COSINE_RBF_H_
#define SRC_RAISED_COSINE_RBF_H_
#include "RBF_Function.h"
#include <cmath>


class Raised_Cosine_RBF : public RBF_Function {
public:
	Raised_Cosine_RBF(MAT A, VEC x, VEC c, FL alpha = 0) :
		RBF_Function(A, x, c, alpha) {}
	Raised_Cosine_RBF(){}

	~Raised_Cosine_RBF() {}

	FL phi_dash(){
		FL pi_by_alpha = M_PI/alpha;
		FL mu = calc_mu();

		FL t1 =  sin(pi_by_alpha * (1 - mu)) * pi_by_alpha * 0.5 * (cos(M_PI * mu) + 1);

		FL t_tmp = cos( M_PI * (1 - mu)/(2*alpha));
		FL t2 = -1 * t_tmp * t_tmp * M_PI *sin(M_PI * mu);

		return t1 + t2;
	}

	FL dphi_dalpha(){
		FL mu = calc_mu();
		FL t1 = M_PI * 0.5 / (alpha * alpha);
		FL t2 = sin(M_PI * (1 - mu)/ alpha );
		FL t3 = cos(M_PI * mu) + 1;
		return t1 * t2 * t3;
	}

	FL get_val(){
		FL mu = calc_mu();
		FL t = cos(M_PI * (1- mu) / (2* alpha));
		FL t2 = cos(M_PI * mu) + 1;
		return t * t * t2;
	}
};



#endif /* SRC_RAISED_COSINE_RBF_H_ */
