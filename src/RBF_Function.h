/*
 * RBF_Function.h
 *
 *  Created on: Jun 6, 2016
 *      Author: chaturvedi
 */

#ifndef SRC_RBF_FUNCTION_H_
#define SRC_RBF_FUNCTION_H_

#include "helpers.h"
using namespace std;

class RBF_Function {
public:
	MAT A;
	VEC x, c;
	FL alpha;
	const FL EPS_MU = 1e-16;
	RBF_Function();

	RBF_Function(MAT A, VEC x, VEC c, FL alpha = 0) {
		this->A = A;
		this->x = x;
		this->c = c;
		this->alpha = alpha;
	}

//	virtual ~RBF_Function() {}

	virtual FL phi_dash() {}
	virtual FL dphi_dalpha(){}
	virtual FL get_val(){}


	FL calc_mu(){
		FL mu;
		mu = sqrt(arma::as_scalar(x_minux_c().t() * this->A * this->A.t() * x_minux_c()));
		if(mu < EPS_MU)
			cout << "mu is zero" <<endl;
		return mu;
	}

	VEC dmu_dc(){
		return -1*(this->A * this->A.t() * x_minux_c())/calc_mu();
	}

	MAT dmu_dA(){
		return A*(x_minux_c()*x_minux_c().t())/calc_mu();
	}

private:
VEC x_minux_c(){
	return this->x - this->c;
}

};

class Shape_Function{

public:
	VEC lambda, x, c;

	Shape_Function();
	Shape_Function(VEC x, VEC c, VEC lambda){
		this->x = x;
		this->c = c;
		this->lambda = lambda;
	}
//	virtual ~Shape_Function() {}

	VEC drho_dc(){
		return -1 * this->lambda;
	}

	VEC drho_dlambda(){
		return (this-> x - this->c);
	}

	FL calc_rho(){
		return as_scalar(this->lambda.t() * (this->x - this->c));
	}

	virtual FL z_dash(){}
	virtual FL get_val(){}
};


#endif /* SRC_RBF_FUNCTION_H_ */
