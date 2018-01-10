#ifndef SRC_LP_simplex_H_
#define SRC_LP_simplex_H_

#include "helpers.h"
// using namespace std;

class LP_simplex{
public:
	LP_simplex() {}
	BOOL B_inv_update_status = false;
	void set_LP(MAT &A, VEC &b, VEC &c){
		this->A = A; // design matrix
		this->b = b; // constraint variable
		this->c = c; // cost variables
	}

	void set_B_inv(MAT &B_inv){
		this->B_inv = B_inv;
		this->B_inv_update_status = true;
	}

	void update_basic_idx(const UVEC &new_basic_idx){
		this->basic_idx = new_basic_idx;
		this->B = A.cols(new_basic_idx);
		this->B_inv_update_status = false;
	}

	MAT get_A(){return A;}
	MAT get_B(){return B;}
	MAT get_B_inv(){return B_inv;}
	VEC get_b(){return b;}
	VEC get_c(){return c;}
	VEC get_x(){return x;}
	UVEC get_basic_idx(){return basic_idx;}


private:
	MAT A,B,B_inv;
	VEC b,c,x;
	UVEC basic_idx;

};

#endif
