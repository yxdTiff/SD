//
//  rc.h
//  sdlb
//
//  Created by Yifan Liu on 6/26/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//

#ifndef sdlb_rc_h
#define sdlb_rc_h

#include "sdglobal.h"
BOOL get_index_number(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s, int obs_idx);
int encode_col(prob_type *p, unsigned long *col, int *cstat, int word_length);
int encode_row(prob_type *p, unsigned long *row, int *rstat, int word_length);
int decode_col(prob_type *p, int *col_num, unsigned long *col, int word_length);
i_type compute_istar_index(sdglobal_type *sd_global, soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                           delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                           double *argmax, BOOL pi_eval, int ictr);
i_type compute_new_istar_index(sdglobal_type *sd_global, soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                               delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                               double *argmax, int ictr);
void new_phi(soln_type *s, prob_type *p, int num_of_phi);
void free_phi(id_type *index);
int get_phi_val(soln_type *s, prob_type *p, id_type *index);
int get_cost_val(sdglobal_type *sd_global, omega_type *omega, num_type *num, id_type *index, double *cost, int *col, int obs_idx);
int adjust_dual_solution(double *Pi, num_type *num, id_type *index);
int calc_theta_from_phi(double *Pi, num_type *num, id_type *index);
int put_phi_into_sigma_delta(sdglobal_type* sd_global, cell_type *c, soln_type *s,
                             lambda_type *lambda, sigma_type *sigma, delta_type *delta,
                             omega_type *omega, num_type *num, sparse_vect *Rbar,
                             sparse_matrix *Tbar,id_type *index);
int adjust_argmax_value(soln_type *s, int obs, sigma_type *sigma, delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,double *arg, id_type *index);
int adjust_alpha_value(soln_type *s, int obs, one_cut *cut ,sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index, int weight);
int adjust_beta_value(soln_type *s, int obs, one_cut *cut ,sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index, int weight);
int adjust_incumbent_height(soln_type *s, int obs, one_cut *cut , double *beta, sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num, id_type *index);
int check_pi_feasibility(sdglobal_type* sd_global, cell_type *c, soln_type *s, prob_type *p,
                         omega_type *omega, num_type *num, ids_type *ids, BOOL new_pi, BOOL new_omega);
BOOL calculate_reduced_cost(sdglobal_type* sd_global, cell_type *c, soln_type *s, prob_type *p,
                            omega_type *omega, num_type *num, ids_type *ids, id_type* index, double *pi);
int reconstruct_pi(sdglobal_type* sd_global, lambda_type *lambda, sigma_type *sigma, delta_type *delta,
                   omega_type *omega, num_type *num, id_type *index, double *pi);
int calc_W_trans_nu(sdglobal_type *sd_global, cell_type *c, soln_type *s, prob_type *p,
                    omega_type *omega, num_type *num, ids_type *ids, id_type *index, double *nu);
#endif
