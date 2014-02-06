//
//  resumeb.h
//  sdlb
//
//  Created by Yifan Liu on 2/4/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//

#ifndef sdlb_resumeb_h
#define sdlb_resumeb_h

#include "sdglobal.h"
int store_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s);
int restore_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s);

#endif
