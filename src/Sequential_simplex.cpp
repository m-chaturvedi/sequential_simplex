/*
 * Sequential_simplex.cpp
 *
 *  Created on: Jan 4, 2018
 *      Author: xiaofengma
 */

#include "Sequential_simplex.h"
#include "LP_simplex.h"

using namespace std;
using namespace arma;

FL Sequential_simplex::get_phi_times_z(){
	FL fx;
	fx = this->phi->get_val()*this->z->get_val();
	return fx;
}

void Sequential_simplex::build_initial_LP(INTEGER center_num){
	INTEGER n = center_num + 1;
	VEC en,e2n;
	MAT En,E2n;
	en = zeros(n,1);
	e2n = zeros(2*n,1);
	En = eye(n,n);
	E2n = eye(2*n,2*n);
	// definition of initial A
	MAT M1,M2,M3,M4,A,A_neg,A_aug,A_initial;
	M1 = En;
	M2 = -1%En;
	M3 = join_horiz(M1,M2);
	M4 = join_horiz(M2,M2);
	A = join_vert(M3,M4);
	A_neg = -1%A;
	A_aug = join_horiz(A,A_neg);
	A_initial = join_horiz(A_aug,E2n);
	// definition of initial b
	VEC b_initial;
	b_initial = zeros(2*n,1);
	// definition of initial c
	VEC c,c1,c2,c3,c4,c_initial;
	c1 = zeros(n,1);
	c2 = ones(n,1);
	c3 = join_vert(c1,c2);
	c4 = -1%c3;
	c = join_vert(c3,c4);
	c_initial = join_vert(c,e2n);

	// definition of initial x
	VEC x;
	x = zeros(6*n,1);
	// definition of initial basic_idx
	UVEC basic_idx;
	basic_idx = regspace<uvec>(0,1,2*n-1);
	// definition of initial B and B_inv
	MAT B,B_inv;
	B = A_initial.cols(0,2*n-1);
	B_inv = inv(B);
	//initial LP setup
	this->old_LP.set_LP(A_initial,b_initial,c_initial);
	this->old_LP.update_basic_idx(basic_idx);
	this->old_LP.set_B_inv(B_inv);
}

void Sequential_simplex::update_LP(const MAT &Phi){
	INTEGER center_num;
	center_num = Phi.n_cols;
	MAT Phi_neg;
	Phi_neg = -1%Phi;
	// build auxi matrix A_aux
	INTEGER n,m;
	n = this->old_LP.A.n_cols;
	m = this->old_LP.A.n_rows;
	MAT A_aux;
	A_aux = zeros(2,n+2);
	A_aux(0,span(0,center_num-1)) = Phi;
	A_aux(0,span(2*center_num,3*center_num-1)) = Phi_neg;
	A_aux(1,span(0,center_num-1)) = Phi_neg;
	A_aux(1,span(2*center_num,3*center_num-1)) = Phi;
	A_aux(0,n) = 1;
	A_aux(1,n+1) = 1;
	// Determine if new constraints are feasible
	FL b1,b2,t1,t2,x1,x2;
	b1 = this->y + this->epsilon;
	b2 = -1*this->y + this->epsilon;
	t1 = as_scalar(A_aux(0,span(0,n-1))*this->old_LP.x);
	t2 = as_scalar(A_aux(1,span(0,n-1))*this->old_LP.x);
	x1 = b1-t1;
	x2 = b2-t2;
	if(x1>-1*EPS && x2>--1*EPS){
		cout<<"new constraint feasible"<<endl;
		this->tmp_LP = this->old_LP;
	}
	else{
		MAT new_A,new_B_inv;
		VEC new_b,new_c,new_x;
		UVEC new_basic_idx;
		//update A
		new_A = zeros(m+2,n+2);
		new_A(span(0,m-1),span(n-1)) = this->old_LP.A;
		new_A.rows(m,m+1) = A_aux;
		//update basic_idx
		new_basic_idx = zeros<uvec>(m+2,1);
		new_basic_idx(span(0,m-1)) = this->old_LP.basic_idx;
		new_basic_idx(m) = n;
		new_basic_idx(m+1) = n+1;
		//update B_inv
		MAT E2,Ax;
		new_B_inv = zeros(m+2,m+2);
		new_B_inv(span(0,m-1),span(0,m-1)) = this->old_LP.B_inv;
		E2 = eye(2,2);
		new_B_inv(span(m,m+1),span(m,m+1)) = E2;
		Ax = -1*A_aux.cols(this->old_LP.basic_idx)*this->old_LP.B_inv;
		new_B_inv(span(m,m+1),span(0,m-1)) = Ax;
		//update b
		new_b = zeros(m+2,1);
		new_b(span(0,m-1)) = this->old_LP.b;
		new_b(m) = b1;
		new_b(m+1) = b2;
		//update c
		new_c = zeros(m+2,1);
		new_c(span(0,m-1)) = this->old_LP.c;
		//update x
		new_x = zeros(m+2,1);
		new_x(span(0,m-1)) = this->old_LP.x;
		new_x(m) = x1;
		new_x(m+1) = x2;
		this->tmp_LP.set_LP(new_A,new_b,new_c);
		this->tmp_LP.update_basic_idx(new_basic_idx);
		this->tmp_LP.set_B_inv(new_B_inv);
	}
}

MAT Sequential_simplex::build_phi_matrix(MAT C,MAT Lambda,MAT A,FL alpha,Sequential_simplex & g){
	INTEGER embed_dim;
	INTEGER center_num;
	VEC current_center;

	embed_dim = C.n_cols;
	center_num = C.n_rows;
	MAT Phi_matrix = arma::zeros(1,1+center_num);
	for(INTEGER i=0;i<center_num;i=i+1){
		current_center = C.row(i).t();
		g.set_c(current_center);
		Phi_matrix(0,i+1) = g.get_phi_times_z();
	}
	Phi_matrix(0,0) = 1;

	return Phi_matrix;
}
