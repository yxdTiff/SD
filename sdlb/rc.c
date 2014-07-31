//
//  rc.c
//  sdlb
//
//  Created by Yifan Liu on 6/26/14.
//  Copyright (c) 2014 Yifan Liu. All rights reserved.
//


#include <time.h>
#include <float.h>
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
#include "resumeb.h"
#include "rc.h"

BOOL get_index_number(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    int i,j;
    int *cstat;
    int *rstat;
    char *basismsg;
    unsigned long *plain, *plain2;
    static int basis_rv_cnt = 0;
    BOOL New_index=TRUE;
    
    /* One extra bit for accomodating index and another for one norm */
    /* Save the 0th location for the norm of the index vector */
    if (!(plain = arr_alloc(s->ids->num_word, unsigned long)))
		err_msg("Allocation", "get_index_number", "plain");
    /* plain2 is for basis status of the slack */
    if (!(plain2 = arr_alloc(s->ids->num_word, unsigned long)))
        err_msg("Allocation", "get_index_number", "plain2");
    
	cstat = (int *) malloc((p->num->sub_cols + 1) * sizeof(int));
    rstat = (int *) malloc((p->num->sub_rows + 1) * sizeof(int));
    
	get_basis(c->subprob, cstat + 1, rstat + 1); /* 2011.10.30 */
    
	if (0)
	{
		/* Write out the solution */
        double *y;
        y = (double *) malloc((p->num->sub_cols + 1) * sizeof(double));
        get_x(c->subprob, y + 1, 0, p->num->sub_cols - 1); /* 2011.10.30 */
		for (j = 1; j <= p->num->sub_cols; j++)
		{
			printf("y[%d] = %17.10g", j, y[j]);
			if (cstat != NULL)
			{
				switch (cstat[j])
				{
                    case AT_LOWER:
                        basismsg = "Nonbasic at lower bound";
                        break;
                    case BASIC:
                        basismsg = "Basic";
                        break;
                    case AT_UPPER:
                        basismsg = "Nonbasic at upper bound";
                        break;
                    case FREE_SUPER:
                        basismsg = "Superbasic, or free variable at zero";
                        break;
                    default:
                        basismsg = "Bad basis status";
                        break;
				}
				printf("  %s   %d", basismsg, cstat[j]);
			}
			printf("\n");
		}
        
	}
    
#if CHECK_PHI
    for (j = 1; j <= p->num->sub_cols; j++) {
        printf("%d", cstat[j]);
    }
    printf("\t %d\n",c->k);
    
    int countc=0,countr=0;
    for (j = 1; j <= p->num->sub_cols; j++) {
        countc = countc + cstat[j];
    }
    printf("countc = %d\t",countc);
    for (j = 1; j <= p->num->sub_rows; j++) {
        countr = countr + rstat[j];
    }
    printf("countr = %d\t sum=%d\n",countr,countc+countr);
#endif

    encode_col(p, plain, cstat, WORD_LENGTH);
    for (j = 1; j < s->ids->num_word; j++) {
        plain[0] += plain[j];
    }
    encode_row(p, plain2, rstat, WORD_LENGTH);
    for (j = 1; j < s->ids->num_word; j++) {
        plain2[0] += plain2[j];
    }
    
    for (i = 0; i < s->ids->cnt; i++) {
        if (equal_ulong_arr(plain, s->ids->index[i]->val, s->ids->num_word) && equal_ulong_arr(plain2, s->ids->index2[i]->val, s->ids->num_word)) {
            New_index = FALSE;
            s->ids->index[i]->freq++;
            s->ids->current_index_idx = i;
        }
    }
    
    if (New_index) {
        s->ids->index[s->ids->cnt] = new_id();
        s->ids->index[s->ids->cnt]->val = plain;
        s->ids->index[s->ids->cnt]->first_c_k = c->k;
        s->ids->index[s->ids->cnt]->freq = 1;
        s->ids->index2[s->ids->cnt] = new_id();
        s->ids->index2[s->ids->cnt]->val = plain2;
        s->ids->index2[s->ids->cnt]->first_c_k = c->k;
        s->ids->index2[s->ids->cnt]->freq = 1;
        s->ids->current_index_idx = s->ids->cnt;
        s->ids->cnt++;
    }
    
#ifdef CHECK_PHI
    printf("%lu\n",s->ids->omega_index[1]);
#endif
    
    for (j = 1; j<=s->ids->num_word; j++) {
        /* Bitwise "B and O" */
        s->rcdata->phi_col[j] = plain[j] & s->ids->omega_index[j];
        /* Bitwise "(not B) and O " */
        s->rcdata->lhs_chl[j] = (~plain[j]) & s->ids->omega_index[j];
    }
    
#ifdef CHECK_PHI
    printf("Here are the Basis columns:\n");
#endif
    decode_col(p, s->rcdata->col_num, plain, WORD_LENGTH);
    /* Let's decode Phi column. Note: phi_col's zero-th location store the actual norm */
#ifdef CHECK_PHI
    printf("Here are the Phi columns:\n");
#endif
    s->rcdata->phi_col[0] = decode_col(p, s->rcdata->phi_col_num, s->rcdata->phi_col, WORD_LENGTH);
    /* Let's decode lhs column. Note: lhs_col's zero-th location store the actual norm */
#ifdef CHECK_PHI
    printf("Here are the Nonbasis random columns:\n");
#endif
    s->rcdata->lhs_chl[0] = decode_col(p, s->rcdata->lhs_col_num, s->rcdata->lhs_chl, WORD_LENGTH);
    if (s->rcdata->phi_col[0]) {
        basis_rv_cnt++;
    }
    
#ifdef CHECK_PHI
    int newhead[2];
    double phi[2];
    int k;
    /* Let's print out all the phi value */
    if (p->num->rv_g) {
        get_basis_head(p->subprob, newhead);
        printf("newhead[%d]=%d\n", 0, newhead[0]);
        printf("newhead[%d]=%d\n", 1, newhead[1]);
        for (i = 0; i < p->num->sub_rows; i++) {
            for (j = 0; j < p->num->rv_g; j++) {
                /* Note the newhead starts with 0 while s->rcdata->phi_col_num starts with 1 */
                if (newhead[i] == s->rcdata->phi_col_num[j]-1) {
                    get_basis_row(p->subprob, i, phi);
                    printf("The following is the phi we want:\n");
                    for(k = 0; k < 2; k++) {
                        printf("phi[%d]:%f\n", k, phi[k]);
                    }
                    printf("\n");
                }
            }
        }
    }
#endif
//    FILE *index_number;
//    index_number = fopen("col.txt", "a");
//    fprintf(index_number, "%f",(double) basis_rv_cnt/c->LP_cnt);
//    fprintf(index_number, "\n");
//    fclose(index_number);

    if (!New_index) {
        mem_free(plain);
    }
    
    mem_free(cstat);
    mem_free(rstat);
    
    return New_index;
}
/* Decode the column and return the total number of 1's */
int decode_col(prob_type *p, int *col_num, unsigned long *col, int word_length)
{
    /* Use unsigned long, so that when shift the bit right, left will be padded with 0's */
    unsigned long mask, temp;
    mask = 1;
    int j, cnt, group, shift;
    /* Let's decode phi_col */
    cnt = 0;
    for (j = 1; j <= p->num->sub_cols; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) col[group] >> shift;
        temp = temp & mask;
        if (temp) {
            col_num[cnt] = j;
#ifdef CHECK_PHI
            printf("column number is %d\n", col_num[cnt]);
#endif
            cnt++;
        }
    }
    if (!cnt) {
#ifdef CHECK_PHI
        printf("None!\n");
#endif
    }
    
    return cnt;
}

int encode_col(prob_type *p, unsigned long *col, int *cstat, int word_length)
{
    int j,group,shift;
    unsigned long temp;
    for (j = 1; j <= p->num->sub_cols; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) cstat[j] << shift;
        col[group] |= temp;
    }
    return 0;
}

int encode_row(prob_type *p, unsigned long *row, int *rstat, int word_length)
{
    int j,group,shift;
    unsigned long temp;
    for (j = 1; j <= p->num->sub_rows; j++) {
        group = j/word_length + 1;
        shift = word_length - j%word_length;
        temp = (unsigned long) rstat[j] << shift;
        row[group] |= temp;
    }
    return 0;
}

/***********************************************************************\
 ** This function loops through all the dual vectors found so far
 ** and returns the index of the one which satisfies the expression:
 **
 **	argmax { Pi x (R - T x X) | all Pi }
 **
 ** where X, R, and T are given.  It is calculated in this form:
 **
 ** 	Pi x Rbar + Pi x Romega + (Pi x Tbar) x X + (Pi x Tomega) x X
 **
 ** Since the Pi's are indexed by the index sets which point to two
 ** different structures (sigma and delta). The following procedure goes
 ** through all the index set and retrieves all sigma and delta needed
 ** for the calculation.
 \***********************************************************************/
i_type compute_istar_index(soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                     delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                     double *argmax, BOOL pi_eval, int ictr)
{
    double arg;
    //  double	argmax;             //modified by Yifan to return argmax value 09/22/2011
    int sig_pi, del_pi, index_idx;
    int c, new_pisz;
    i_type ans;
    
#ifdef LOOP
    printf("Inside compute_istar_index\n");
#endif
    //printf("\n***IN COMPUTE ISTAR ...\n");
    if (pi_eval == TRUE)
    {
        //new_pisz =sd_global->config.PI_EVAL_START;
        new_pisz = ictr / 10 + 1;
    }
    else
        new_pisz = 0;
    
    ictr -= new_pisz; /*evaluate the pi's generated in the first 90% iterations */
    
    *argmax = -DBL_MAX;
    
    /*added by Yifan to enable parallel process*/
    // #pragma omp for private(sig_pi,del_pi, arg)
    for (index_idx = 0 ; index_idx < s->ids->cnt; index_idx++) {
        sig_pi = s->ids->sig_idx[index_idx];
        if (sigma->ck[sig_pi] <= ictr) {
            /* Find the row in delta corresponding to this row in sigma */
            del_pi = sigma->lamb[sig_pi];
            
            /* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
            arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R
            - Pi_Tbar_X[sig_pi];
            
            /* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
            
            for (c = 1; c <= num->rv_cols; c++)
                arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];
            
#ifdef LOOP
            print_sigma(sigma, num, sig_pi);
            print_delta(delta, num, del_pi, obs);
            printf("\nResulting arg=%f", arg);
#endif
            
            if (arg > (*argmax))
            {
                *argmax = arg;
                ans.sigma = sig_pi;
                ans.delta = del_pi;
                //printf("argmax:%f and istar(%d,%d)\n", arg, ans.sigma, ans.delta);
            }
        }
    }
    
#ifdef LOOP
    printf("Exiting compute_istar_index\n");
#endif
    
    return ans;
}

i_type compute_new_istar_index(soln_type *s, int obs, one_cut *cut, sigma_type *sigma,
                         delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
                         double *argmax, int ictr)
{
    double arg;
    int sig_pi, del_pi, index_idx;
    int c, new_pisz;
    i_type ans;
    ans.sigma = 0;
    ans.delta = 0;
    
#ifdef LOOP
    printf("Inside compute_new_istar_index\n");
#endif
    
    new_pisz = ictr / 10 + 1;
    ictr -= new_pisz; /*evaluate the pi's generated in the last 10% iterations */
    
    *argmax = -DBL_MAX;
    
    for (index_idx = 0; index_idx < s->ids->cnt; index_idx++) {
        sig_pi = s->ids->sig_idx[index_idx];
        if (sigma->ck[sig_pi] > ictr) {
            /* Find the row in delta corresponding to this row in sigma */
            del_pi = sigma->lamb[sig_pi];
            
            /* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
            arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R
            - Pi_Tbar_X[sig_pi];
            
            /* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
            for (c = 1; c <= num->rv_cols; c++)
                arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];
            
#ifdef LOOP
            print_sigma(sigma, num, sig_pi);
            print_delta(delta, num, del_pi, obs);
            printf("\nResulting arg=%f", arg);
#endif
            
            if (arg > (*argmax))
            {
                *argmax = arg;
                ans.sigma = sig_pi;
                ans.delta = del_pi;
                
            }
        }
    }
    
#ifdef LOOP
    printf("Exiting compute_new_istar_index\n");
#endif
    return ans;
}

