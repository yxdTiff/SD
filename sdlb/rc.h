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
int decode_col(prob_type *p, int *col_num, unsigned long *col, int word_length);
#endif
