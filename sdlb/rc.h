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
BOOL get_index_number(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s);
int encode_col(prob_type *p, unsigned long *col, int *cstat, int word_length);
int encode_row(prob_type *p, unsigned long *row, int *rstat, int word_length);
int decode_col(prob_type *p, int *col_num, unsigned long *col, int word_length);
i_type compute_istar_index(soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                           delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                           double *argmax, BOOL pi_eval, int ictr);
i_type compute_new_istar_index(soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                               delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                               double *argmax, int ictr);
#endif
