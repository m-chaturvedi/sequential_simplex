/*
 * sequential_simplex.h
 *
 *  Created on: Jan 3, 2018
 *      Author: xiaofengma
 */

#include "RBF_Function.h"
#include "LP_simplex.h"

class Sequential_simplex{
public:
	LP_simplex old_LP,new_LP,tmp_LP;
	void set_functions(RBF_Function &phi,Shape_Function &z){
		this->phi = &phi;
		this->z = &z;
	}
	void set_x(const VEC &p_x){this->x = p_x; this->phi->x = p_x; this->z->x = p_x;}
	void set_y(const FL &y){this->y = y;}
	void set_c(const VEC &p_c){this->c = p_c; this->phi->c = p_c; this->z->c = p_c;}
	void set_A(const MAT &p_A){this->phi->A = p_A;}
	void set_lambda(const VEC &p_lambda){this->lambda = p_lambda; this->z->lambda = p_lambda;}
	void set_alpha(const FL &p_alpha){this->alpha = p_alpha; this->phi->alpha = p_alpha;}
	//TODO:may need to change epsilon while running
	void set_epsilon(const FL &p_epsilon){this->epsilon = p_epsilon;}

	FL get_phi_times_z();
	void build_initial_LP(INTEGER center_num);
	void update_LP(const MAT &Phi);

	static MAT build_phi_matrix(MAT C,MAT Lambda,MAT A,FL alpha, Sequential_simplex &g);

private:
	RBF_Function *phi;
	Shape_Function *z;
	VEC x;
	VEC c;
	VEC lambda;
	FL alpha;
	FL epsilon;
	FL y;

};
