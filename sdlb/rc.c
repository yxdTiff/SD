//
//  rc.c
//  sdlb
//
//  Created by Yifan Liu on 6/26/14.
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
#include "resumeb.h"
#include "rc.h"

int get_index_number(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
    FILE *index_number;
    int j;
    int *cstat;
    int *rstat;
    int sum = 0;
    double *y;
    char *basismsg;
	cstat = (int *) malloc((p->num->sub_cols + 1) * sizeof(int));
    rstat = (int *) malloc((p->num->sub_rows + 1) * sizeof(int));
	y = (double *) malloc((p->num->sub_cols + 1) * sizeof(double));
    
	get_basis(c->subprob, cstat + 1, rstat + 1); /* 2011.10.30 */
	get_x(c->subprob, y + 1, 0, p->num->sub_cols - 1); /* 2011.10.30 */
    
	if (0)
	{
		/* Write out the solution */
        
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
    index_number = fopen("indexNumber.txt", "a");
    for (j = 1; j <= p->num->sub_cols; j++) {
        fprintf(index_number, "%d", cstat[j]);
        sum = sum + cstat[j];
    }
    for (j = 1; j <= p->num->sub_rows; j++) {
        fprintf(index_number, "%d", rstat[j]);
        sum = sum + rstat[j];
    }
    printf("%d\n",sum);
    fprintf(index_number, "\n");
    fclose(index_number);
    
    return 0;
}