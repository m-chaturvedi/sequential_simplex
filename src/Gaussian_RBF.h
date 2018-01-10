/*
 * Gaussian_RBF.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: chaturvedi
 */
#ifndef SRC_GAUSSIAN_RBF_H_
#define SRC_GAUSSIAN_RBF_H_

#include "RBF_Function.h"

class Gaussian_RBF : public RBF_Function{
public:
	Gaussian_RBF(MAT A, VEC x, VEC c, FL alpha = 0) :
		RBF_Function(A, x, c, alpha) {}
	Gaussian_RBF(){}

	~Gaussian_RBF() {}

	FL phi_dash(){
		return (-2*calc_mu()/(alpha*alpha))*exp(-(calc_mu() * calc_mu())/(alpha * alpha));
	}

	FL dphi_dalpha(){
		return (2*calc_mu()*calc_mu()/ (alpha * alpha *alpha))* exp(-(calc_mu()*calc_mu())/(alpha * alpha));
	}

	FL get_val(){
		return exp(-calc_mu()*calc_mu()/ (alpha * alpha));
	}
};
#endif

