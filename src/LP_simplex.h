#ifndef SRC_LP_simplex_H_
#define SRC_LP_simplex_H_

#include "helper.h"
using namespace std;

class LP_simplex{
public:
	void set_A(MAT &A, VEC &b, VEC &c){
		this->A = A; // design matrix
		this->b = b; // constraint variable
		this->c = c; // cost variables
	}

	MAT B;
	VEC b,c,x;
	UVEC basic_idx;

	void update_lp(const UVEC &new_basic_idx){
		this->basic_idx = new_basic_idx;
		this->B = A.cols(new_basic_idx);
	}


	BOOL feasibility_test(){

	}

private:
	MAT A;
};

#endif
