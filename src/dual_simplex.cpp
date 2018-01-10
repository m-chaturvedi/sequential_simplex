/*
 * dual_simplex.cpp
 *
 *  Created on: Oct 22, 2017
 *      Author: xiaofengma
 */
#include "dual_simplex.h"
//#include "LP_simplex.h"

using namespace std;
using namespace arma;

dual_simplex::dual_simplex(LP_simplex &my_LP){
	this->my_LP = &my_LP;
	this->A = my_LP.get_A();
	this->b = my_LP.get_b();
	this->c = my_LP.get_c();
	this->old_basic_idx = my_LP.get_basic_idx();
	this->B_inv = my_LP.get_B_inv();
	this->optimal_flag = false;
	this->l = -1;
	this->j = -1;
}

/*compute_c_bar compute the reduced cost c_bar*/
void dual_simplex::compute_c_bar(){
	// c_bar = c'-c_B'*B_inv*A
	VEC c_B = this->c.rows(this->old_basic_idx);
	this->c_bar = this->c.t() - c_B.t() * this->B_inv * this->A;
	assert_msg(all(this->c_bar>-EPS), " c_bar should be non negative i.e. solution should be optimal");
}

/*compute_l() find the leaving index and update the optimality flag*/
void dual_simplex::compute_l(){
	VEC x_B = this->B_inv * this->b;
	if(all(x_B > -EPS)){
		this->optimal_flag = true;
	}
	else{
		//find the first l such that X_{B(l)} < 0
		UVEC negative_idx = find(x_B < -EPS);
		this->l = negative_idx(0);
	}
}

/*compute_i() find the entering index j by looking at the negative entries in B(l)-th row*/
void dual_simplex::compute_j(){
	ROWVEC v = B.row(this->l) * this->A;
	//UVEC negatvie_idx = find(v < -EPS);
	if(all(v < -EPS)){
		assert_msg(false,"optimal dual cost is infinity");
	}
	else{
		UVEC negative_idx = find(v < -EPS);
		VEC v_neg = v(negative_idx);
		ROWVEC c_bar_neg = this->c_bar(negative_idx);
		VEC ratio = c_bar_neg.t() / abs(v_neg);
		this->j = ratio.index_min(); // not sure if index_min() will always return the smallest index?
	}
}

/*compute_new_basic_idx() update the new basic index given i and j*/
void dual_simplex::compute_new_basic_idx(){
	this->new_basic_idx = this->old_basic_idx;
	this->new_basic_idx(this->l) = j;
}

/*compute_B_inv() compute the inverse of B based on the old B_inv*/
void dual_simplex::compute_B_inv(){
	uword m = this->B_inv.n_rows;
	UVEC u = this->B_inv * this->A.col(this->j);
	MAT Q = eye<mat>(m,m);
	UVEC Q_l = -u/u(this->l);
	Q_l(this->l) = -Q_l(l)/u(l);
	Q.row(this->l) = Q_l;
	this->B_inv = Q*this->B_inv;
}

/*run_dual_simplex() iterates the dual simplex until optimal solution or max iteration reached*/
void dual_simplex::run_dual_simplex(){
	INTEGER count = 0;
	while(count < MAX_ITERS || this->optimal_flag == false){
		this->compute_c_bar();
		this->compute_l();
		if(this->optimal_flag == true){
			this->my_LP->update_basic_idx(this->new_basic_idx);
			this->my_LP->set_B_inv(this->B_inv);
			break;
		}
		else{
		this->compute_j();
		this->compute_new_basic_idx();
		this->compute_B_inv();
		}
		count += count;
		if(count==MAX_ITERS-1){
			cout<<"maximum number of iterations reached! Discard this point!"<<endl;
		}
	}
}
