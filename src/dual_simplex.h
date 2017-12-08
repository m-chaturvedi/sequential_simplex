/*
 * dual_simplex.h
 *
 *  Created on: Oct 21, 2017
 *      Author: xiaofengma
 */

#ifndef SRC_DUAL_SIMPLEX_H_
#define SRC_DUAL_SIMPLEX_H_

#include "helper.h"
#include "LP_simplex.h"

class dual_simplex{
	// take in a LP with basic non-feasible solution solve dual for feasibility
public:
	dual_simplex(LP_simplex &my_LP);
	dual_simplex(){}
	void compute_l();
	void compute_j();
	void compute_c_bar();
	void compute_B_inv();
	void compute_new_basic_idx();
	void run_dual_simplex();
	void update_LP();
	bool get_optimal_flag(){
		return this->optimal_flag;
	};

private:
	LP_simplex * my_LP;
	ROWVEC c_bar;
	VEC c;
	VEC b;
	MAT B_inv;
	MAT A,B;
	arma::uword l; // choose l such that X_{B(l)} < 0
	arma::uword j;// for all v_i < 0, compute \overline{c}_i / \abs{v_i}, j is the index corresponding to the smallest ratio
	UVEC new_basic_idx;
	UVEC old_basic_idx;
	bool optimal_flag;
};



#endif /* SRC_DUAL_SIMPLEX_H_ */
