/*
 * omega.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef OMEGA_H_
#define OMEGA_H_
#include "sdglobal.h"

BOOL equal_obs(int *a, int *b, int len);
BOOL valid_omega_idx(omega_type *omega, int idx);
int generate_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num,
		BOOL *new_omeg, sd_long *RUN_SEED);
int get_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num, BOOL *new_omeg);
int next_omega_idx(omega_type *omega);
omega_type *new_omega(int num_iter, int num_rv, coord_type *coord);
ids_type *new_ids(int num_iter, int sub_col);
id_type *new_id(void);
rc_type *new_rcdata(int sub_rows, int rv_g, int num_word);
int get_omega_index(prob_type *p, soln_type *s);
void free_omega(omega_type *omega);
void free_ids(ids_type *ids);
void free_id(id_type *index);
void free_rcdata(rc_type *rc_data);
void get_R_T_G_omega(sdglobal_type* sd_global, omega_type *omega, int obs_idx);
void init_R_T_G_omega(sparse_vect *Romega, sparse_matrix *Tomega, sparse_vect *Gomega,
		omega_type *omega, num_type *num);
void print_omega(omega_type *omega, num_type *num, int idx);

#endif /* OMEGA_H_ */
