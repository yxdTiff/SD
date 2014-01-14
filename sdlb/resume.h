//
//  resume.h
//  sdlb
//
//  Created by Yifan Liu on 1/12/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//

#ifndef sdlb_resume_h
#define sdlb_resume_h
#include "sdglobal.h"
int store_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s);
int restore_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s);
#endif
