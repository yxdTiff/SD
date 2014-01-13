//
//  resume.c
//  sdlb
//
//  Created by Yifan Liu on 1/12/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//

#include <time.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "theta.h"
#include "optimal.h"
#include "omega.h"
#include "log.h"
#include "master.h"
#include "cuts.h"
#include "rvgen.h"
#include "sdconstants.h"
#include "sdglobal.h"
#include "resume.h"

int store_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int cnt, idx, obs;
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    
    /* 1. Start storing omega_type */
    printf("\n");
    printf("OMEGA\n");
    printf("cnt: %d\n",s->omega->cnt);
    printf("next: %d\n",s->omega->next);
    printf("most: %d\n",s->omega->most);
    printf("last: %d\n",s->omega->last);
    printf("k: %d\n",s->omega->k);
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        printf("row[%d]: %d\n", cnt, s->omega->row[cnt]);
    }
    for (cnt = 0; cnt <= p->num->rv; cnt++) {
        printf("col[%d]: %d\n", cnt, s->omega->col[cnt]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->cipher; cnt++) {
            printf("idx[%d][%d]: %d\n", idx, cnt, s->omega->idx[idx][cnt]);
        }
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        printf("weight[%d]: %d\n", idx, s->omega->weight[idx]);
    }
    for (idx = 0; idx < s->omega->cnt; idx++) {
        printf("filter[%d]: %d\n", idx, s->omega->filter[idx]);
    }
    for (cnt = 1; cnt <= sd_global->omegas.num_omega; cnt++) {
        printf("RT[%d]: %f\n", cnt, s->omega->RT[cnt]);
    }
    
    /* 2. Start storing lambda_type */
    printf("LAMBDA\n");
    printf("cnt: %d\n",c->lambda->cnt);
    for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
        printf("row[%d]: %d\n", cnt, c->lambda->row[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (cnt = 0; cnt <= p->num->rv_rows; cnt++) {
            printf("val[%d][%d]: %f\n", idx, cnt, c->lambda->val[idx][cnt]);
        }
    }
    
    /* 3. Start storing sigma_type */
    printf("SIGMA\n");
    printf("cnt: %d\n",c->sigma->cnt);
    for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
        printf("col[%d]: %d\n", cnt, c->sigma->col[cnt]);
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        printf("val[%d].R: %f\n", idx, c->sigma->val[idx].R);
        for (cnt = 0; cnt <= p->num->nz_cols; cnt++) {
            printf("val[%d].T[%d]: %f\n", idx, cnt, c->sigma->val[idx].T[cnt]);
        }
    }
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        printf("ck[%d]: %d\n", idx, c->sigma->ck[idx]);
    }
    
    /* 4. Start storing theta_type */
    printf("TEHTA\n");
    printf("cnt: %d\n",c->theta->cnt);
    for (idx = 0; idx < c->theta->cnt; idx++) {
        printf("last[%d]: %d\n", idx, c->theta->last[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        printf("k[%d]: %f\n", idx, c->theta->k[idx]);
    }
    for (idx = 0; idx < c->theta->cnt; idx++) {
        printf("p[%d]: %f\n", idx, c->theta->p[idx]);
    }
    
    /* 5. Start storing delta_type */
    printf("DELTA\n");
    for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
        printf("col[%d]: %d\n", cnt, s->delta->col[cnt]);
    }
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                printf("val[%d][%d].R: %f\n", idx, obs, s->delta->val[idx][obs].R);
                for (cnt = 0; cnt <= p->num->rv_cols; cnt++) {
                    printf("val[%d][%d].T[%d]: %f\n", idx, obs, cnt, s->delta->val[idx][obs].T[cnt]);
                }
            }
        }
    }

    /* 6. Start storing cuts */
    printf("CUTS\n");
    
    /*  */
    
    return 0;
}