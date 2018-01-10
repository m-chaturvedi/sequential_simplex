/*
 * test.cpp
 *
 *  Created on: Jan 3, 2018
 *      Author: xiaofengma
 */

#include "helpers.h"
#include "RBF_Function.h"
#include "Raised_Cosine_RBF.h"

using namespace arma;

int main(){
	MAT A = zeros(3,3);
	VEC x = zeros(3,1);
	VEC c = zeros(3,1);
	FL alpha = 0;
	Raised_Cosine_RBF a(A,x,c,alpha);
	return 0;
}

