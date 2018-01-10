/*
 * Shape_Function.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: chaturvedi
 */

#ifndef SRC_ARC_TAN_SHAPE_FUNCTION_H_
#define SRC_ARC_TAN_SHAPE_FUNCTION_H_
#include <cmath>
#include "RBF_Function.h"

class Arc_Tan : public Shape_Function {
public:
	Arc_Tan(VEC x, VEC c, VEC lambda) : Shape_Function(x, c, lambda) {}
	~Arc_Tan() {}

	FL get_val(){
		return atan(this->calc_rho())/M_PI + 0.5;
	}

	FL z_dash(){
		return 1/(M_PI * (1 + this->calc_rho() * this->calc_rho()));
	}
};

#endif


