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
    printf("cnt: %d\n", c->cuts->cnt);
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        printf("val[%d]->omega_cnt: %d\n",idx, c->cuts->val[idx]->omega_cnt);
        for (cnt = 0; cnt < c->cuts->val[idx]->omega_cnt; cnt++) {
            printf("val[%d]->istar[%d]: %d\n", idx, cnt, c->cuts->val[idx]->istar[cnt]);
        }
        printf("val[%d]->slack_cnt: %d\n",idx, c->cuts->val[idx]->slack_cnt);
        printf("val[%d]->cell_num: %d\n",idx, c->cuts->val[idx]->cell_num);
        printf("val[%d]->cut_obs: %d\n",idx, c->cuts->val[idx]->cut_obs);
        printf("val[%d]->row_num: %d\n",idx, c->cuts->val[idx]->row_num);
        printf("val[%d]->alpha: %f\n",idx, c->cuts->val[idx]->alpha);
        printf("val[%d]->alpha_incumb: %f\n",idx, c->cuts->val[idx]->alpha_incumb);
        for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
            printf("val[%d]->beta[%d]: %f\n", idx, cnt, c->cuts->val[idx]->beta[cnt]);
        }
        printf("val[%d]->subfeaflag: %u\n",idx, c->cuts->val[idx]->subfeaflag);
        printf("val[%d]->is_incumbent: %u\n",idx, c->cuts->val[idx]->is_incumbent);
        if (c->cuts->val[idx]->is_incumbent) {
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                printf("val[%d]->subobj_omega[%d]: %f\n", idx, obs, c->cuts->val[idx]->subobj_omega[obs]);
            }
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                printf("val[%d]->subobj_omega[%d]: %d\n", idx, obs, c->cuts->val[idx]->subobj_freq[obs]);
            }
        }
        
    }
    
    /* 6. Start storing run_time */
    printf("TIME\n");
    printf("total_time: %f\n", s->run_time->total_time);
    printf("iteration_time: %f\n", s->run_time->iteration_time);
    printf("soln_master_iter: %f\n", s->run_time->soln_master_iter);
    printf("soln_subprob_iter: %f\n", s->run_time->soln_subprob_iter);
    printf("full_test_iter: %f\n", s->run_time->full_test_iter);
    printf("argmax_iter: %f\n", s->run_time->argmax_iter);
    printf("iteration_accum: %f\n", s->run_time->iteration_accum);
    printf("soln_master_accum: %f\n", s->run_time->soln_master_accum);
    printf("soln_subprob_accum: %f\n", s->run_time->soln_subprob_accum);
    printf("full_test_accum: %f\n", s->run_time->full_test_accum);
    printf("argmax_accum: %f\n", s->run_time->argmax_accum);



    /* 8. Start storing the remaining of cell structure */
    printf("CELL\n");
    printf("id_num: %d\n", c->id_num);
    printf("num_members: %d\n", c->num_members);
    for (idx = 0; idx < c->num_members; idx++) {
        printf("members[%d]: %d\n", idx, c->members[idx]);
    }
    printf("quad_scalar: %f\n", c->quad_scalar);
    printf("LP_cnt: %d\n", c->LP_cnt);
    printf("LP_test: %d\n", c->LP_test);
    printf("N: %d\n", c->N);
    printf("P: %f\n", c->P);
    printf("k: %d\n", c->k);
    printf("opt_mode: %u\n", c->opt_mode);
    printf("incumb_infea: %u\n", c->incumb_infea);
    printf("fea_count: %d\n", c->fea_count);
    
    /* 9. Start storing the remaining of soln structure */
    printf("SOLN\n");
    for (idx = 0; idx <= p->num->sub_rows; idx++) {
        printf("Pi[%d]: %f\n", idx, s->Pi[idx]);
    }
    printf("subobj_est: %f\n", s->subobj_est);
    for (idx = 0; idx <= p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt; idx++) {
        printf("Master_pi[%d]: %f\n", idx, s->Master_pi[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("Master_dj[%d]: %f\n", idx, s->Master_dj[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("candid_x[%d]: %f\n", idx, s->candid_x[idx]);
    }
    printf("candid_est: %f\n", s->candid_est);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_x[%d]: %f\n", idx, s->incumb_x[idx]);
    }
    printf("incumb_k: %d\n", s->incumb_k);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_d[%d]: %f\n", idx, s->incumb_d[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_avg[%d]: %f\n", idx, s->incumb_avg[idx]);
    }
    printf("alpha: %f\n", s->alpha);
    for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
        printf("beta[%d]: %f\n", cnt, s->beta[cnt]);
    }
    printf("incumb_est: %f\n", s->incumb_est);
    printf("opt_value: %f\n", s->opt_value);
    printf("norm_d_k_1: %f\n", s->norm_d_k_1);
    printf("norm_d_k: %f\n", s->norm_d_k);
    printf("incumb_stdev: %f\n", s->incumb_stdev);
    printf("incumb_cut: %d\n", s->incumb_cut);
    printf("last_update: %d\n", s->last_update);
    printf("gamma: %f\n", s->gamma);
    printf("optimality_flag: %u\n", s->optimality_flag);
    printf("smpl_ever_flag: %u\n", s->smpl_ever_flag);
    printf("smpl_test_flag: %u\n", s->smpl_test_flag);
    printf("incumbent_change: %u\n", s->incumbent_change);
    printf("dual_statble_flag: %u\n", *(s->dual_statble_flag));
    printf("full_test_error: %f\n", s->full_test_error);
    printf("passed: %d\n", s->passed);
    printf("sub_lb_checker: %f\n", s->sub_lb_checker);
    printf("max_ratio: %f\n", s->max_ratio);
    printf("min_ratio: %f\n", s->min_ratio);
    for (idx = 0; idx < sd_global->config.MAX_SCAN_LEN; idx++) {
        printf("pi_ratio[%d]: %f\n", idx, s->pi_ratio[idx]);
    }

    /* 9. Start storing the last seed used in getting observation */
    printf("config.RUN_SEED: %lld\n", sd_global->config.RUN_SEED);
    
    return 0;
}

int restore_sd_data(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
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
    printf("cnt: %d\n", c->cuts->cnt);
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        printf("val[%d]->omega_cnt: %d\n",idx, c->cuts->val[idx]->omega_cnt);
        for (cnt = 0; cnt < c->cuts->val[idx]->omega_cnt; cnt++) {
            printf("val[%d]->istar[%d]: %d\n", idx, cnt, c->cuts->val[idx]->istar[cnt]);
        }
        printf("val[%d]->slack_cnt: %d\n",idx, c->cuts->val[idx]->slack_cnt);
        printf("val[%d]->cell_num: %d\n",idx, c->cuts->val[idx]->cell_num);
        printf("val[%d]->cut_obs: %d\n",idx, c->cuts->val[idx]->cut_obs);
        printf("val[%d]->row_num: %d\n",idx, c->cuts->val[idx]->row_num);
        printf("val[%d]->alpha: %f\n",idx, c->cuts->val[idx]->alpha);
        printf("val[%d]->alpha_incumb: %f\n",idx, c->cuts->val[idx]->alpha_incumb);
        for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
            printf("val[%d]->beta[%d]: %f\n", idx, cnt, c->cuts->val[idx]->beta[cnt]);
        }
        printf("val[%d]->subfeaflag: %u\n",idx, c->cuts->val[idx]->subfeaflag);
        printf("val[%d]->is_incumbent: %u\n",idx, c->cuts->val[idx]->is_incumbent);
        if (c->cuts->val[idx]->is_incumbent) {
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                printf("val[%d]->subobj_omega[%d]: %f\n", idx, obs, c->cuts->val[idx]->subobj_omega[obs]);
            }
            for (obs = 0; obs < c->cuts->val[idx]->omega_cnt; obs++) {
                printf("val[%d]->subobj_omega[%d]: %d\n", idx, obs, c->cuts->val[idx]->subobj_freq[obs]);
            }
        }
        
    }
    
    /* 6. Start storing run_time */
    printf("TIME\n");
    printf("total_time: %f\n", s->run_time->total_time);
    printf("iteration_time: %f\n", s->run_time->iteration_time);
    printf("soln_master_iter: %f\n", s->run_time->soln_master_iter);
    printf("soln_subprob_iter: %f\n", s->run_time->soln_subprob_iter);
    printf("full_test_iter: %f\n", s->run_time->full_test_iter);
    printf("argmax_iter: %f\n", s->run_time->argmax_iter);
    printf("iteration_accum: %f\n", s->run_time->iteration_accum);
    printf("soln_master_accum: %f\n", s->run_time->soln_master_accum);
    printf("soln_subprob_accum: %f\n", s->run_time->soln_subprob_accum);
    printf("full_test_accum: %f\n", s->run_time->full_test_accum);
    printf("argmax_accum: %f\n", s->run_time->argmax_accum);
    
    
    
    /* 8. Start storing the remaining of cell structure */
    printf("CELL\n");
    printf("id_num: %d\n", c->id_num);
    printf("num_members: %d\n", c->num_members);
    for (idx = 0; idx < c->num_members; idx++) {
        printf("members[%d]: %d\n", idx, c->members[idx]);
    }
    printf("quad_scalar: %f\n", c->quad_scalar);
    printf("LP_cnt: %d\n", c->LP_cnt);
    printf("LP_test: %d\n", c->LP_test);
    printf("N: %d\n", c->N);
    printf("P: %f\n", c->P);
    printf("k: %d\n", c->k);
    printf("opt_mode: %u\n", c->opt_mode);
    printf("incumb_infea: %u\n", c->incumb_infea);
    printf("fea_count: %d\n", c->fea_count);
    
    /* 9. Start storing the remaining of soln structure */
    printf("SOLN\n");
    for (idx = 0; idx <= p->num->sub_rows; idx++) {
        printf("Pi[%d]: %f\n", idx, s->Pi[idx]);
    }
    printf("subobj_est: %f\n", s->subobj_est);
    for (idx = 0; idx <= p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt; idx++) {
        printf("Master_pi[%d]: %f\n", idx, s->Master_pi[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("Master_dj[%d]: %f\n", idx, s->Master_dj[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("candid_x[%d]: %f\n", idx, s->candid_x[idx]);
    }
    printf("candid_est: %f\n", s->candid_est);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_x[%d]: %f\n", idx, s->incumb_x[idx]);
    }
    printf("incumb_k: %d\n", s->incumb_k);
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_d[%d]: %f\n", idx, s->incumb_d[idx]);
    }
    for (idx = 0; idx <= p->num->mast_cols; idx++) {
        printf("incumb_avg[%d]: %f\n", idx, s->incumb_avg[idx]);
    }
    printf("alpha: %f\n", s->alpha);
    for (cnt = 0; cnt <= p->num->mast_cols; cnt++) {
        printf("beta[%d]: %f\n", cnt, s->beta[cnt]);
    }
    printf("incumb_est: %f\n", s->incumb_est);
    printf("opt_value: %f\n", s->opt_value);
    printf("norm_d_k_1: %f\n", s->norm_d_k_1);
    printf("norm_d_k: %f\n", s->norm_d_k);
    printf("incumb_stdev: %f\n", s->incumb_stdev);
    printf("incumb_cut: %d\n", s->incumb_cut);
    printf("last_update: %d\n", s->last_update);
    printf("gamma: %f\n", s->gamma);
    printf("optimality_flag: %u\n", s->optimality_flag);
    printf("smpl_ever_flag: %u\n", s->smpl_ever_flag);
    printf("smpl_test_flag: %u\n", s->smpl_test_flag);
    printf("incumbent_change: %u\n", s->incumbent_change);
    printf("dual_statble_flag: %u\n", *(s->dual_statble_flag));
    printf("full_test_error: %f\n", s->full_test_error);
    printf("passed: %d\n", s->passed);
    printf("sub_lb_checker: %f\n", s->sub_lb_checker);
    printf("max_ratio: %f\n", s->max_ratio);
    printf("min_ratio: %f\n", s->min_ratio);
    for (idx = 0; idx < sd_global->config.MAX_SCAN_LEN; idx++) {
        printf("pi_ratio[%d]: %f\n", idx, s->pi_ratio[idx]);
    }
    
    /* 9. Start storing the last seed used in getting observation */
    printf("config.RUN_SEED: %lld\n", sd_global->config.RUN_SEED);
    
    return 0;
}