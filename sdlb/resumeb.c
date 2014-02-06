//
//  resumeb.c
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

int store_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int idx, obs;
    FILE *rep_data;
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen("rep_data.txt", "w");
    /* 1. Start storing omega_type */
    fwrite(&(s->omega->cnt), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->next), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->most), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->last), sizeof(int), 1, rep_data);
    fwrite(&(s->omega->k), sizeof(int), 1, rep_data);

    fwrite(s->omega->row, sizeof(int), p->num->rv+1, rep_data);
    
    
    fwrite(s->omega->col, sizeof(int), p->num->rv+1, rep_data);
    
    for (idx = 0; idx < s->omega->cnt; idx++) {

        fwrite(s->omega->idx[idx], sizeof(int), p->num->cipher+1, rep_data);
        
    }

    fwrite(s->omega->weight, sizeof(int), s->omega->cnt, rep_data);
    

    fwrite(s->omega->filter, sizeof(int), s->omega->cnt, rep_data);

    
    fwrite(s->omega->RT+1, sizeof(double), sd_global->omegas.num_omega, rep_data);
    
    
    /* 2. Start storing lambda_type */
    fwrite(&(c->lambda->cnt), sizeof(int), 1, rep_data);

    fwrite(c->lambda->row, sizeof(int), p->num->rv_rows+1, rep_data);
    
    for (idx = 0; idx < c->lambda->cnt; idx++) {

        fwrite(c->lambda->val[idx], sizeof(double), p->num->rv_rows+1, rep_data);
        
    }
    
    /* 3. Start storing sigma_type */
    fwrite(&(c->sigma->cnt), sizeof(int), 1, rep_data);

    fwrite(c->sigma->col, sizeof(int), p->num->nz_cols+1, rep_data);
    
    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fwrite(&(c->sigma->val[idx].R), sizeof(double), 1, rep_data);

        fwrite(c->sigma->val[idx].T, sizeof(double), p->num->nz_cols+1, rep_data);
        
    }

    fwrite(c->sigma->ck, sizeof(int), c->sigma->cnt, rep_data);
    

    fwrite(c->sigma->lamb, sizeof(int), c->sigma->cnt, rep_data);
    
    
    /* 4. Start storing theta_type */
    fwrite(&(c->theta->cnt), sizeof(int), 1, rep_data);

    fwrite(c->theta->last, sizeof(int), c->theta->cnt, rep_data);
    

    fwrite(c->theta->k, sizeof(double), c->theta->cnt, rep_data);
    

    fwrite(c->theta->p, sizeof(double), c->theta->cnt, rep_data);
    
    
    /* 5. Start storing delta_type */

    fwrite(s->delta->col, sizeof(int), p->num->rv_cols+1, rep_data);
    
    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                fwrite(&(s->delta->val[idx][obs].R), sizeof(double), 1, rep_data);

                fwrite(s->delta->val[idx][obs].T, sizeof(double), p->num->rv_cols+1, rep_data);
                
            }
        }
    }
    
    /* 6. Start storing cuts */
    fwrite(&(c->cuts->cnt), sizeof(int), 1, rep_data);
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        fwrite(&(c->cuts->val[idx]->omega_cnt), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->cut_obs), sizeof(int), 1, rep_data);

        fwrite(c->cuts->val[idx]->istar, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
        
        fwrite(&(c->cuts->val[idx]->slack_cnt), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->cell_num), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->row_num), sizeof(int), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->alpha), sizeof(double), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->alpha_incumb), sizeof(double), 1, rep_data);

        fwrite(c->cuts->val[idx]->beta, sizeof(double), p->num->mast_cols+1, rep_data);
        
        fwrite(&(c->cuts->val[idx]->subfeaflag), sizeof(BOOL), 1, rep_data);
        fwrite(&(c->cuts->val[idx]->is_incumbent), sizeof(BOOL), 1, rep_data);
        if (c->cuts->val[idx]->is_incumbent) {

            fwrite(c->cuts->val[idx]->subobj_omega, sizeof(double), c->cuts->val[idx]->omega_cnt, rep_data);
            

            fwrite(c->cuts->val[idx]->subobj_freq, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
            
        }
        
    }
    
    /* 7. Start storing run_time */
    fwrite(&(s->run_time->total_time), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->iteration_time), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_master_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_subprob_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->full_test_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->argmax_iter), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->iteration_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_master_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->soln_subprob_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->full_test_accum), sizeof(double), 1, rep_data);
    fwrite(&(s->run_time->argmax_accum), sizeof(double), 1, rep_data);
    
    
    
    /* 8. Start storing the remaining of cell structure */
    fwrite(&(c->id_num), sizeof(int), 1, rep_data);
    fwrite(&(c->num_members), sizeof(int), 1, rep_data);

    fwrite(c->members, sizeof(int), c->num_members, rep_data);
    
    fwrite(&(c->quad_scalar), sizeof(double), 1, rep_data);
    fwrite(&(c->LP_cnt), sizeof(int), 1, rep_data);
    fwrite(&(c->LP_test), sizeof(int), 1, rep_data);
    fwrite(&(c->N), sizeof(int), 1, rep_data);
    fwrite(&(c->P), sizeof(double), 1, rep_data);
    fwrite(&(c->k), sizeof(int), 1, rep_data);
    fwrite(&(c->opt_mode), sizeof(BOOL), 1, rep_data);
    fwrite(&(c->incumb_infea), sizeof(BOOL), 1, rep_data);
    fwrite(&(c->fea_count), sizeof(int), 1, rep_data);
    
    /* 9. Start storing the remaining of soln structure */

    fwrite(s->Pi, sizeof(double), p->num->sub_rows+1, rep_data);
    
    fwrite(&(s->subobj_est), sizeof(double), 1, rep_data);

    fwrite(s->Master_pi, sizeof(double), p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt + 1, rep_data);
    

    fwrite(s->Master_dj, sizeof(double), p->num->mast_cols + 1, rep_data);
    

    fwrite(s->candid_x, sizeof(double), p->num->mast_cols + 1, rep_data);

    fwrite(&(s->candid_est), sizeof(double), 1, rep_data);

    fwrite(s->incumb_x, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fwrite(&(s->incumb_k), sizeof(int), 1, rep_data);

    fwrite(s->incumb_d, sizeof(double), p->num->mast_cols + 1, rep_data);
    

    fwrite(s->incumb_avg, sizeof(double), p->num->mast_cols + 1, rep_data);
 
    fwrite(&(s->alpha), sizeof(double), 1, rep_data);

    fwrite(s->beta, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fwrite(&(s->incumb_est), sizeof(double), 1, rep_data);
    fwrite(&(s->opt_value), sizeof(double), 1, rep_data);
    fwrite(&(s->norm_d_k_1), sizeof(double), 1, rep_data);
    fwrite(&(s->norm_d_k), sizeof(double), 1, rep_data);
    fwrite(&(s->incumb_stdev), sizeof(double), 1, rep_data);
    fwrite(&(s->incumb_cut), sizeof(int), 1, rep_data);
    fwrite(&(s->last_update), sizeof(int), 1, rep_data);
    fwrite(&(s->gamma), sizeof(double), 1, rep_data);
    fwrite(&(s->optimality_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->smpl_ever_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->smpl_test_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->incumbent_change), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->dual_statble_flag), sizeof(BOOL), 1, rep_data);
    fwrite(&(s->full_test_error), sizeof(double), 1, rep_data);
    fwrite(&(s->passed), sizeof(int), 1, rep_data);
    fwrite(&(s->sub_lb_checker), sizeof(double), 1, rep_data);
    fwrite(&(s->max_ratio), sizeof(double), 1, rep_data);
    fwrite(&(s->min_ratio), sizeof(double), 1, rep_data);

    fwrite(s->pi_ratio, sizeof(double), sd_global->config.MAX_SCAN_LEN, rep_data);
    
    
    /* 10. Start storing the last seed used in getting observation */
    fwrite(&(sd_global->config.RUN_SEED), sizeof(sd_long), 1, rep_data);
    
    fclose(rep_data);
    
    return 0;
}

int restore_sd_data_b(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int cnt, idx, obs;
    FILE *rep_data;
    /* Store all data structure necessary for resuming SD */
    /* modified by Yifan 2014.01.12 */
    rep_data = fopen("rep_data.txt", "r");
    
    /* 1. Start restoring omega_type */
    fread(&(s->omega->cnt), sizeof(int), 1, rep_data);
    for (idx = 0; idx < s->omega->cnt; idx++) {
        s->omega->idx[idx] = arr_alloc(p->num->cipher+1, int);
    }
    
    fread(&(s->omega->next), sizeof(int), 1, rep_data);
    fread(&(s->omega->most), sizeof(int), 1, rep_data);
    fread(&(s->omega->last), sizeof(int), 1, rep_data);
    fread(&(s->omega->k), sizeof(int), 1, rep_data);
    
    fread(s->omega->row, sizeof(int), p->num->rv+1, rep_data);


    fread(s->omega->col, sizeof(int), p->num->rv+1, rep_data);
    
    for (idx = 0; idx < s->omega->cnt; idx++) {
        
        fread(s->omega->idx[idx], sizeof(int), p->num->cipher+1, rep_data);
        
    }
    
    fread(s->omega->weight, sizeof(int), s->omega->cnt, rep_data);
    
    
    fread(s->omega->filter, sizeof(int), s->omega->cnt, rep_data);
    
    
    fread(s->omega->RT+1, sizeof(double), sd_global->omegas.num_omega, rep_data);
    
    
    /* 2. Start restoring lambda_type */
    fread(&(c->lambda->cnt), sizeof(int), 1, rep_data);
    
    /* Dynamically allocate memory space for lambda */
    for (cnt = 0; cnt < c->lambda->cnt; cnt++)
		if (!(c->lambda->val[cnt] = arr_alloc(p->num->rv_rows+1, double)))
			err_msg("Allocation", "restore_lambda", "lambda->val[cnt]");
    

    fread(c->lambda->row, sizeof(int), p->num->rv_rows+1, rep_data);

    for (idx = 0; idx < c->lambda->cnt; idx++) {
        
        fread(c->lambda->val[idx], sizeof(double), p->num->rv_rows+1, rep_data);
        
    }
    
    /* 3. Start restoring sigma_type */
    fread(&(c->sigma->cnt), sizeof(int), 1, rep_data);
    
    /* Dynamically allocate memory space for sigma */
    for (cnt = 0; cnt < c->sigma->cnt && cnt < sd_global->config.MAX_ITER; cnt++)
		if (!(c->sigma->val[cnt].T = arr_alloc(p->num->nz_cols+1, double)))
			err_msg("Allocation", "resotre_sigma", "sigma->val[cnt]");
    
    
    fread(c->sigma->col, sizeof(int), p->num->nz_cols+1, rep_data);

    for (idx = 0; idx < c->sigma->cnt; idx++) {
        fread(&(c->sigma->val[idx].R), sizeof(double), 1, rep_data);
        
        fread(c->sigma->val[idx].T, sizeof(double), p->num->nz_cols+1, rep_data);
        
    }
    
    fread(c->sigma->ck, sizeof(int), c->sigma->cnt, rep_data);

    
    fread(c->sigma->lamb, sizeof(int), c->sigma->cnt, rep_data);

    
    /* 4. Start restoring theta_type */
    fread(&(c->theta->cnt), sizeof(int), 1, rep_data);
    
    fread(c->theta->last, sizeof(int), c->theta->cnt, rep_data);

    
    fread(c->theta->k, sizeof(double), c->theta->cnt, rep_data);

    
    fread(c->theta->p, sizeof(double), c->theta->cnt, rep_data);

    
    
    /* 5. Start restoring delta_type */
    
    /* Dynamically allocate memory space for delta */
    /* Double check this part since there might be an error here!!!! */
    for (cnt = 0; cnt < c->lambda->cnt; cnt++){
        if (!(s->delta->val[cnt] = arr_alloc(p->num->iter, pi_R_T_type)))
            err_msg("Allocation", "restore_delta", "delta->val");
        for (obs=0; obs<s->omega->most; obs++) {
            if (!(s->delta->val[cnt][obs].T = arr_alloc(p->num->rv_cols + 1, double)))
                err_msg("Allocation", "restore_delta", "delta->val[cnt]");
        }
    }
    
    
    
    fread(s->delta->col, sizeof(int), p->num->rv_cols+1, rep_data);

    for (idx = 0; idx < c->lambda->cnt; idx++) {
        for (obs = 0; obs < s->omega->most; obs++) {
            if (valid_omega_idx(s->omega, obs)) {
                fread(&(s->delta->val[idx][obs].R), sizeof(double), 1, rep_data);
                
                fread(s->delta->val[idx][obs].T, sizeof(double), p->num->rv_cols+1, rep_data);
                
            }
        }
    }
    
    /* 6. Start restoring cuts */
    fread(&(c->cuts->cnt), sizeof(int), 1, rep_data);
    
    /* Dynamically allocate memory space for cuts */
    
    int omega_cnt, cut_obs;
    for (idx = 0; idx < c->cuts->cnt; idx++) {
        fread(&(omega_cnt), sizeof(int), 1, rep_data);
        fread(&(cut_obs), sizeof(int), 1, rep_data);
        c->cuts->val[idx] = new_cut(p->num->mast_cols, omega_cnt, cut_obs);
        c->cuts->val[idx]->omega_cnt = omega_cnt;
        c->cuts->val[idx]->cut_obs = cut_obs;
        
        fread(c->cuts->val[idx]->istar, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
        
        fread(&(c->cuts->val[idx]->slack_cnt), sizeof(int), 1, rep_data);
        fread(&(c->cuts->val[idx]->cell_num), sizeof(int), 1, rep_data);
        fread(&(c->cuts->val[idx]->row_num), sizeof(int), 1, rep_data);
        fread(&(c->cuts->val[idx]->alpha), sizeof(double), 1, rep_data);
        fread(&(c->cuts->val[idx]->alpha_incumb), sizeof(double), 1, rep_data);
        
        fread(c->cuts->val[idx]->beta, sizeof(double), p->num->mast_cols+1, rep_data);
        
        fread(&(c->cuts->val[idx]->subfeaflag), sizeof(BOOL), 1, rep_data);
        fread(&(c->cuts->val[idx]->is_incumbent), sizeof(BOOL), 1, rep_data);
        if (c->cuts->val[idx]->is_incumbent) {
            if (!(c->cuts->val[idx]->subobj_omega = arr_alloc( c->cuts->val[idx]->omega_cnt+1, double)))
                err_msg("Allocation", "SD_cut", "subobj_omega");
            if (!(c->cuts->val[idx]->subobj_freq = arr_alloc( c->cuts->val[idx]->omega_cnt+1, int)))
                err_msg("Allocation", "SD_cut", "subobj_freq");
            
            fread(c->cuts->val[idx]->subobj_omega, sizeof(double), c->cuts->val[idx]->omega_cnt, rep_data);
            
            
            fread(c->cuts->val[idx]->subobj_freq, sizeof(int), c->cuts->val[idx]->omega_cnt, rep_data);
            
        }
        
    }
    
    /* 7. Start restoring run_time */
    fread(&(s->run_time->total_time), sizeof(double), 1, rep_data);
    fread(&(s->run_time->iteration_time), sizeof(double), 1, rep_data);
    fread(&(s->run_time->soln_master_iter), sizeof(double), 1, rep_data);
    fread(&(s->run_time->soln_subprob_iter), sizeof(double), 1, rep_data);
    fread(&(s->run_time->full_test_iter), sizeof(double), 1, rep_data);
    fread(&(s->run_time->argmax_iter), sizeof(double), 1, rep_data);
    fread(&(s->run_time->iteration_accum), sizeof(double), 1, rep_data);
    fread(&(s->run_time->soln_master_accum), sizeof(double), 1, rep_data);
    fread(&(s->run_time->soln_subprob_accum), sizeof(double), 1, rep_data);
    fread(&(s->run_time->full_test_accum), sizeof(double), 1, rep_data);
    fread(&(s->run_time->argmax_accum), sizeof(double), 1, rep_data);
    
    
    
    /* 8. Start restoring the remaining of cell structure */
    fread(&(c->id_num), sizeof(int), 1, rep_data);
    fread(&(c->num_members), sizeof(int), 1, rep_data);
    
    fread(c->members, sizeof(int), c->num_members, rep_data);
    
    fread(&(c->quad_scalar), sizeof(double), 1, rep_data);
    fread(&(c->LP_cnt), sizeof(int), 1, rep_data);
    fread(&(c->LP_test), sizeof(int), 1, rep_data);
    fread(&(c->N), sizeof(int), 1, rep_data);
    fread(&(c->P), sizeof(double), 1, rep_data);
    fread(&(c->k), sizeof(int), 1, rep_data);
    fread(&(c->opt_mode), sizeof(BOOL), 1, rep_data);
    fread(&(c->incumb_infea), sizeof(BOOL), 1, rep_data);
    fread(&(c->fea_count), sizeof(int), 1, rep_data);
    
    /* 9. Start restoring the remaining of soln structure */
    
    fread(s->Pi, sizeof(double), p->num->sub_rows+1, rep_data);
    
    fread(&(s->subobj_est), sizeof(double), 1, rep_data);
    /* Update the dual size before restoring added by Yifan 2014.01.28*/
    update_dual_size(c, s, p);
    
    fread(s->Master_pi, sizeof(double), p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt + 1, rep_data);
    

    
    fread(s->Master_dj, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    
    fread(s->candid_x, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fread(&(s->candid_est), sizeof(double), 1, rep_data);
    
    fread(s->incumb_x, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fread(&(s->incumb_k), sizeof(int), 1, rep_data);
    
    fread(s->incumb_d, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    
    fread(s->incumb_avg, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fread(&(s->alpha), sizeof(double), 1, rep_data);
    
    fread(s->beta, sizeof(double), p->num->mast_cols + 1, rep_data);
    
    fread(&(s->incumb_est), sizeof(double), 1, rep_data);
    fread(&(s->opt_value), sizeof(double), 1, rep_data);
    fread(&(s->norm_d_k_1), sizeof(double), 1, rep_data);
    fread(&(s->norm_d_k), sizeof(double), 1, rep_data);
    fread(&(s->incumb_stdev), sizeof(double), 1, rep_data);
    fread(&(s->incumb_cut), sizeof(int), 1, rep_data);
    fread(&(s->last_update), sizeof(int), 1, rep_data);
    fread(&(s->gamma), sizeof(double), 1, rep_data);
    fread(&(s->optimality_flag), sizeof(BOOL), 1, rep_data);
    fread(&(s->smpl_ever_flag), sizeof(BOOL), 1, rep_data);
    fread(&(s->smpl_test_flag), sizeof(BOOL), 1, rep_data);
    fread(&(s->incumbent_change), sizeof(BOOL), 1, rep_data);
    fread((s->dual_statble_flag), sizeof(BOOL), 1, rep_data);
    /* make sure dual_stable_flag is turned to NO */
    *s->dual_statble_flag = 0;
    fread(&(s->full_test_error), sizeof(double), 1, rep_data);
    fread(&(s->passed), sizeof(int), 1, rep_data);
    fread(&(s->sub_lb_checker), sizeof(double), 1, rep_data);
    fread(&(s->max_ratio), sizeof(double), 1, rep_data);
    fread(&(s->min_ratio), sizeof(double), 1, rep_data);
    
    fread(s->pi_ratio, sizeof(double), sd_global->config.MAX_SCAN_LEN, rep_data);
    
    
    /* 10. Start restoring the last seed used in getting observation */
    fread(&(sd_global->config.RUN_SEED), sizeof(sd_long), 1, rep_data);
    
    fclose(rep_data);
    
    return 0;
}