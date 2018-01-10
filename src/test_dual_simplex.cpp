/*
 * test_dual_simplex.cpp
 *
 *  Created on: Jan 9, 2018
 *      Author: xiaofengma
 */
#include <iostream>
#include "helpers.h"
#include "RBF_Function.h"
#include "Arc_Tan_Shape_Function.h"
#include "Raised_Cosine_RBF.h"
#include "LP_simplex.h"
#include "Dual_simplex.h"

using namespace std;
using namespace arma;

int main(){
	LP_simplex test_LP;

	MAT A,B_inv;
	VEC b;
	VEC c;
	UVEC basic_idx;

	A<<-1<<-2<<1<<0<<endr
	 <<-1<<0<<0<<1<<endr;

	B_inv<<1<<0<<endr
		 <<0<<1<<endr;

	b<<-2<<endr
	 <<-1<<endr;

	c<<1<<endr
	 <<1<<endr
	 <<0<<endr
	 <<0<<endr;

	basic_idx<<3<<endr
			 <<4<<endr;

	test_LP.set_LP(A,b,c);
	test_LP.update_basic_idx(basic_idx);
	test_LP.set_B_inv(B_inv);

	dual_simplex g(test_LP);

	return 0;
}



