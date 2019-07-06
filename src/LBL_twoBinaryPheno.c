#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <time.h>
#include <omp.h>
#include <assert.h>

void update_beta1_bi(int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int which_beta1,int N, double lambda, double exp_u[N], double sigma_sq_u, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat, int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1]);
void update_beta2_bi(int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int which_beta2, int N, double lambda, double exp_u[N], double sigma_sq_u, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat,  int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1]);
double update_lambda_bi(double *beta1, double *beta2, double a, double b, int x_length);
void update_u(int N, double *u, double exp_u[N], double sigma_sq_u, int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat, int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1]);
double update_sigma_sq_u(double sigma_sq_u, double *u, int N);
double update_D_bi(double *freq1, double *freq2, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, int *y1_new, int *y2_new, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map);
void update_freq_00(double *freq_00, double *freq_10, double *freq_01, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map);
void update_freq_10(double *freq_00, double *freq_10, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map);
void update_freq_01(double *freq_00, double *freq_01, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map);

double calc_a_bi(double *freq, int *per_xz, double D);
double sum_bi(double *x, int n);
void dirichlet_bi(double *param, int dim, double *gen_sample);
double find_min_bi(double *arr, int n);
double gen_double_exp_bi(double mean, double SD);

void mcmc_twoBinaryPheno(int *haplotype_map, int *tot_hap, int *y1, int *y2, int *N, int *n00, int *n10, int *n01, int *num_haplo_id, int *x_length, double *freq_00, double *freq_10, double *freq_01, double *D, double *beta1, double *beta2, double *a, double *b, int *unique_map, int *tot_uniq_mat, double *lambda, double *u, double *sigma_sq_u, int *NUM_IT, int *BURN_IN, double *beta1_out, double *beta2_out, double *lambda_out, double *freq_00_out, double *freq_10_out, double *freq_01_out, double *D_out, double *sigma_sq_u_out)
{
	/*  sum_indepmary of notations
	y1[i]: responce(case/control) of disease 1 for each subject (i<*N)
	y2[i]: responce(case/control) of disease 2 for each subject (i<*N)

	*N: num of subjects
	*tot_hap: nrow of the design matrix
	num_haplo_id[i]: num of rows for each subject (i<*N)

	*x_length: length of design matrix (length of beta vector=*x_length+1)

	haplo_map[i][j][k]: map of haplotypes (i<*N, j<num_haplo_id[i], k<2)
	uniq_map[i][j]: map of haplotypes
	*tot_uniq_mat: num of unique rows of haplo.map
	*/
	int i, j, k, l, m, n, first, second, which_beta, it = 0, it1 = 0, it2 = 0/*,it3 = 0*/;
	int ***haplo_map, h_mat[*tot_hap][2], **uniq_map;
  int y1_new[*N], y2_new[*N];
	double freq_1[*x_length + 1], freq_2[*x_length + 1], temp1, temp2;
	double exp_Xbeta1_map[*x_length + 1][*x_length + 1], exp_Xbeta2_map[*x_length + 1][*x_length + 1], a1_map[*x_length + 1][*x_length + 1], a2_map[*x_length + 1][*x_length + 1], exp_u[*N];
  double **theta_over_deno_map;

	/*y1_new y2_new*/
	//y1_new = (int *)calloc(*N, sizeof(int));
	//y2_new = (int *)calloc(*N, sizeof(int));

	/* separating haplotype map vector from R as h_mat matrix */
	l = 0;
	for (j = 0; j<2; ++j)
	{
		for (i = 0; i<*tot_hap; ++i)
		{
			h_mat[i][j] = haplotype_map[l];
			++l;
		}
	}

	/* separating h_mat as haplo_map */
	haplo_map = calloc(*N, sizeof(int*));
	for (i = 0; i<*N; ++i)
	{
		haplo_map[i] = calloc(num_haplo_id[i], sizeof(int*));
	}

	for (i = 0; i<*N; ++i)
	{
		for (j = 0; j < num_haplo_id[i]; ++j)
		{
			haplo_map[i][j] = calloc(2, sizeof(int));
		}
	}

	l = 0;
	for (i = 0; i<*N; ++i)
	{
		y1_new[i] = y1[l];
		y2_new[i] = y2[l];
		for (j = 0; j < num_haplo_id[i]; ++j)
		{
			m = 0;
			for (k = 0; k < 2; ++k)
			{
				haplo_map[i][j][k] = h_mat[l][m];
				++m;
			}
			++l;
		}
	}

	/* separating unique_map vector from R as uniq_map matrix */

	uniq_map = calloc(*tot_uniq_mat, sizeof(int*));
	for (i = 0; i<*tot_uniq_mat; ++i)
	{
		uniq_map[i] = calloc(2, sizeof(int));
	}

	l = 0;
	for (i = 0; i<*tot_uniq_mat; ++i)
	{
		for (j = 0; j<2; ++j)
		{
			uniq_map[i][j] = unique_map[l];
			++l;
		}
	}
   /* allocate space for theta_over_deno_map */
  theta_over_deno_map = calloc(*N, sizeof(double*));
	for (i = 0; i<*N; ++i)
	{
		theta_over_deno_map[i] = calloc(num_haplo_id[i], sizeof(double));
	}


  /* Calculate freq_1, freq_2, Xbeta1_map, Xbeta2_map--used in the first iteration to update beta1*/
    for (i = 0; i < *x_length + 1; ++i)
		{
			freq_1[i] = (freq_00[i] * (*n00) + freq_01[i] * (*n01)) / (*n00 + *n01);
			freq_2[i] = (freq_00[i] * (*n00) + freq_10[i] * (*n10)) / (*n00 + *n10);
		}

		for (j = 0; j < *tot_uniq_mat; ++j)
		{
      first = uniq_map[j][0] - 1;
      second = uniq_map[j][1] - 1;
      a1_map[first][second] = calc_a_bi(freq_1, uniq_map[j], *D);
			a2_map[first][second] = calc_a_bi(freq_2, uniq_map[j], *D);
	    temp1 = beta1[0];
	    temp2 = beta2[0];
      if (first < *x_length  && second < *x_length)
		    {
		    	temp1 += beta1[first + 1] + beta1[second + 1];
			    temp2 += beta2[first + 1] + beta2[second + 1];
		    }
		    else if (first == *x_length && second < *x_length)
		    {
    			temp1 += beta1[second + 1];
		    	temp2 += beta2[second + 1];
		    }
		    else if (first < *x_length && second == *x_length)
		    {
		    	temp1 += beta1[first + 1];
		    	temp2 += beta2[first + 1];
	    	}
		    exp_Xbeta1_map[first][second] = exp(temp1);
		    exp_Xbeta2_map[first][second] = exp(temp2);
    	}
    for (i = 0; i<*N; ++i)
	{
		exp_u[i] = exp(u[i]);
	}
	/*---------------------start MCMC here------------------------------*/
	for (n = 0; n<*NUM_IT; ++n)
	{

    /* update beta1 parameters */
		for (i = 0; i<*x_length + 1; ++i)
		{
			which_beta = i;
      update_beta1_bi(haplo_map, uniq_map, beta1, beta2, which_beta, *N, *lambda, exp_u, *sigma_sq_u, *x_length, y1_new, y2_new, *tot_uniq_mat, num_haplo_id, a1_map, a2_map, exp_Xbeta1_map, exp_Xbeta2_map);
      update_beta2_bi(haplo_map, uniq_map, beta1, beta2, which_beta, *N, *lambda, exp_u, *sigma_sq_u, *x_length, y1_new, y2_new, *tot_uniq_mat, num_haplo_id, a1_map, a2_map, exp_Xbeta2_map, exp_Xbeta1_map);
		}

		/* update lambda */
		*lambda = update_lambda_bi(beta1, beta2, *a, *b, *x_length);

		/* Update u */
    update_u(*N, u, exp_u, *sigma_sq_u, haplo_map, uniq_map, beta1, beta2, *x_length, y1_new, y2_new, *tot_uniq_mat, num_haplo_id, a1_map, a2_map, exp_Xbeta1_map, exp_Xbeta2_map);

		/*Update sigma_sq_u */
		*sigma_sq_u = 	update_sigma_sq_u(*sigma_sq_u, u, *N);

		/* update D */
		/* calculate theta_over_deno_map at the same time -- used in updatings of freq's and D */
    *D = update_D_bi(freq_1, freq_2, *D, *x_length, *N, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, y1_new, y2_new, exp_Xbeta1_map, exp_Xbeta2_map, exp_u, a1_map, a2_map, theta_over_deno_map);

		/* Updating freq */
    update_freq_00(freq_00, freq_10, freq_01, freq_1, freq_2, *D, *x_length, *N, *n00, *n10, *n01, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, exp_Xbeta1_map, exp_Xbeta2_map, exp_u, a1_map, a2_map, theta_over_deno_map);
    update_freq_10(freq_00, freq_10, freq_1, freq_2, *D, *x_length, *N, *n00, *n10, *n01, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, exp_Xbeta2_map, exp_u, a1_map, a2_map, theta_over_deno_map);
    update_freq_01(freq_00, freq_01, freq_1, freq_2, *D, *x_length, *N, *n00, *n10, *n01, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, exp_Xbeta1_map, exp_u, a1_map, a2_map, theta_over_deno_map);

		/* outputs */
		if (n >= *BURN_IN)
		{
			for (i = 0; i<*x_length + 1; ++i)
			{
				beta1_out[it] = beta1[i];
				beta2_out[it] = beta2[i];
				++it;
			}

			lambda_out[it2] = *lambda;

			for (i = 0; i<*x_length + 1; ++i)
			{
				freq_00_out[it1] = freq_00[i];
				freq_10_out[it1] = freq_10[i];
				freq_01_out[it1] = freq_01[i];
				++it1;
			}

			D_out[it2] = *D;
			sigma_sq_u_out[it2] = *sigma_sq_u;
			++it2;
		}
	}
	/*----------------------finish MCMC-----------------------------------------*/

	for (i = 0; i < *N; ++i)
  {
    for (j = 0; j < num_haplo_id[i]; ++j)
    {
      free(haplo_map[i][j]);
    }
    free(haplo_map[i]);
    free(theta_over_deno_map[i]);
  }
	for (i = 0; i < *tot_uniq_mat; ++i) free(uniq_map[i]);
	free(uniq_map);
  free(haplo_map);
  free(theta_over_deno_map);

}

/* update beta1 with fixed Xbeta1_map_old, Xbeta2_map */
void update_beta1_bi(int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int which_beta1,int N, double lambda, double exp_u[N], double sigma_sq_u, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat, int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1])
{
	double beta1_new, g_old = 0, g_new = 0, f_old, f_new, beta1_new_vec[x_length + 1], SD, accept_prob, term[tot_uniq_mat], term_new[tot_uniq_mat];
	double inprod,inprod_new, temp1_new, theta1, theta1_new, theta2, dat, dat_new;
  double exp_Xbeta1_map_new[x_length + 1][x_length + 1];
	int i, j, r, x, first, second;

	beta1_new = gen_double_exp_bi(beta1[which_beta1], sqrt(fabs(beta1[which_beta1])));

	for (i = 0; i<x_length + 1; ++i)
		beta1_new_vec[i] = beta1[i];
	beta1_new_vec[which_beta1] = beta1_new;

	/*calculate Xbeta1_map_new*/
  #pragma omp parallel for private(first, second, temp1_new)
	for (j = 0; j < tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		temp1_new = beta1_new_vec[0];
		if (first < x_length  && second < x_length)
		{
			temp1_new += beta1_new_vec[first + 1] + beta1_new_vec[second + 1];
		}
		else if (first == x_length && second < x_length)
		{
			temp1_new += beta1_new_vec[second + 1];
		}
		else if (first < x_length && second == x_length)
		{
			temp1_new += beta1_new_vec[first + 1];
		}
		exp_Xbeta1_map_new[first][second] = exp(temp1_new); /* update Xbeta1_map_new here, already have Xbeta1_map_old and Xbeta2_map */
	}

	/*g_old and g_new*/
	g_old = -lambda*fabs(beta1[which_beta1]);
	g_new = -lambda*fabs(beta1_new);

  #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first,second,inprod,inprod_new,term,term_new,theta1,theta1_new,theta2,dat,dat_new)
	for (i = 0; i < N; ++i)
	{

		for (j = 0; j < tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;

			term[j] = a1_map[first][second] * exp_Xbeta1_map[first][second]*exp_u[i];
			term_new[j] = a1_map[first][second] * exp_Xbeta1_map_new[first][second]*exp_u[i];
		}
		g_old += -0.5*log(1 + sum_bi(term, tot_uniq_mat));
		g_new += -0.5*log(1 + sum_bi(term_new, tot_uniq_mat));

    /*Y1=0 & Y2=0, i from 1 to N1 = n_00*/
		if ((y1_new[i] == 0) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second]*exp_u[i];
				theta1_new = exp_Xbeta1_map_new[first][second]*exp_u[i];
				theta2 = exp_Xbeta2_map[first][second]*exp_u[i];
				dat = a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2 + theta1_new*theta2);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);

				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=1 & Y2=0, i from n_00+1 to N2(N2=n_00+n_10)*/
		if ((y1_new[i] == 1) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second]*exp_u[i];
				theta1_new = exp_Xbeta1_map_new[first][second]*exp_u[i];
				theta2 = exp_Xbeta2_map[first][second]*exp_u[i];
				dat = pow(theta1, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta1_new, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2 + theta1_new*theta2);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);

				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);


		}

		/*Y1=0 & Y2=1, i from N2+1 to N3(N3=n_00+n_10+n_01)*/

		if ((y1_new[i] == 0) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second]*exp_u[i];
				theta1_new = exp_Xbeta1_map_new[first][second]*exp_u[i];
				theta2 = exp_Xbeta2_map[first][second]*exp_u[i];
				dat = pow(theta2, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta2, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2 + theta1_new*theta2);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);

				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);

		}

		/*Y1=1 & Y2=1, i from N3+1 to N*/
		if ((y1_new[i] == 1) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second]*exp_u[i];
				theta1_new = exp_Xbeta1_map_new[first][second]*exp_u[i];
				theta2 = exp_Xbeta2_map[first][second]*exp_u[i];
				dat = pow((theta1*theta2), 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow((theta1_new*theta2), 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2 + theta1_new*theta2);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);

				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}
	} /*finish g_old and g_new*/


	  /*f_old and f_new*/
	SD = sqrt(fabs(beta1_new));
	f_old = exp(-sqrt(2)*fabs(beta1[which_beta1] - beta1_new) / SD) / (sqrt(2)*SD);
	SD = sqrt(fabs(beta1[which_beta1]));
	f_new = exp(-sqrt(2)*fabs(beta1[which_beta1] - beta1_new) / SD) / (sqrt(2)*SD);
	accept_prob = exp(g_new - g_old)*(f_old / f_new);

	if (accept_prob > 1)
	{
		beta1[which_beta1] = beta1_new;
		for (j = 0; j < tot_uniq_mat; ++j)
	    {
		    first = uniq_map[j][0] - 1;
		    second = uniq_map[j][1] - 1;
		    exp_Xbeta1_map[first][second] = exp_Xbeta1_map_new[first][second];
	    }
	}
	else
	{
		GetRNGstate();
		x = rbinom(1, accept_prob);
		PutRNGstate();
		if (x == 1)
		{
			beta1[which_beta1] = beta1_new;
      for (j = 0; j < tot_uniq_mat; ++j)
	    {
		    first = uniq_map[j][0] - 1;
		    second = uniq_map[j][1] - 1;
		    exp_Xbeta1_map[first][second] = exp_Xbeta1_map_new[first][second];
	    }
		}
	}
}

void update_beta2_bi(int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int which_beta2, int N, double lambda, double exp_u[N], double sigma_sq_u, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat,  int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1])
{

	double beta2_new, g_old = 0, g_new = 0, f_old, f_new, beta2_new_vec[x_length + 1], SD, accept_prob, term[tot_uniq_mat], term_new[tot_uniq_mat];
	double inprod, inprod_new, temp2_new, theta1, theta2_new, theta2, dat, dat_new;
	double exp_Xbeta2_map_new[x_length + 1][x_length + 1];
	int i, j, r, x, first, second;

	beta2_new = gen_double_exp_bi(beta2[which_beta2], sqrt(fabs(beta2[which_beta2])));


	for (i = 0; i<x_length + 1; ++i)
		beta2_new_vec[i] = beta2[i];
	beta2_new_vec[which_beta2] = beta2_new;

	/*calculate Xbeta1_map, Xbeta1_new_map, Xbeta2_map*/
  #pragma omp parallel for private(first, second, temp2_new)
	for (j = 0; j < tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		temp2_new = beta2_new_vec[0];
		if (first < x_length  && second < x_length)
		{
			temp2_new += beta2_new_vec[first + 1] + beta2_new_vec[second + 1];
		}
		else if (first == x_length && second < x_length)
		{
			temp2_new += beta2_new_vec[second + 1];
		}
		else if (first < x_length && second == x_length)
		{
			temp2_new += beta2_new_vec[first + 1];
		}
		exp_Xbeta2_map_new[first][second] = exp(temp2_new); /* update Xbeta2_map_new here, already have Xbeta2_map_old and Xbeta1_map */
	}

	/*g_old and g_new*/
	g_old = -lambda*fabs(beta2[which_beta2]);
	g_new = -lambda*fabs(beta2_new);

  #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first,second,inprod,inprod_new,term,term_new,theta1,theta2,theta2_new,dat,dat_new)
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
			term[j] = a2_map[first][second] * exp_Xbeta2_map[first][second] * exp_u[i];
			term_new[j] = a2_map[first][second] * exp_Xbeta2_map_new[first][second] * exp_u[i];
		}
		g_old += -0.5*log(1 + sum_bi(term, tot_uniq_mat));
		g_new += -0.5*log(1 + sum_bi(term_new, tot_uniq_mat));
		/*Y1=0 & Y2=0, i from 1 to n_00*/
		if ((y1_new[i] == 0) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta2_new = exp_Xbeta2_map_new[first][second] * exp_u[i];
				dat = a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2_new + theta1*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=1 & Y2=0, i from n_00+1 to N2(N2=n_00+n_10)*/
		if ((y1_new[i] == 1) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta2_new = exp_Xbeta2_map_new[first][second] * exp_u[i];
				dat = pow(theta1, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta1, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2_new + theta1*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=0 & Y2=1, i from N2+1 to N3(N3=n_00+n_10+n_01)*/
		if ((y1_new[i] == 0) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta2_new = exp_Xbeta2_map_new[first][second] * exp_u[i];
				dat = pow(theta2, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta2_new, 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2_new + theta1*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=1 & Y2=1, i from N3+1 to N*/
		if ((y1_new[i] == 1) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta2_new = exp_Xbeta2_map_new[first][second] * exp_u[i];
				dat = pow((theta1*theta2), 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow((theta1*theta2_new), 2)*a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2_new + theta1*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}
	} /*finish g_old and g_new*/

	  /*f_old and f_new*/
	SD = sqrt(fabs(beta2_new));
	f_old = exp(-sqrt(2)*fabs(beta2[which_beta2] - beta2_new) / SD) / (sqrt(2)*SD);
	SD = sqrt(fabs(beta2[which_beta2]));
	f_new = exp(-sqrt(2)*fabs(beta2[which_beta2] - beta2_new) / SD) / (sqrt(2)*SD);
	accept_prob = exp(g_new - g_old)*(f_old / f_new);
	if (accept_prob > 1)
	{
		beta2[which_beta2] = beta2_new;
		for (j = 0; j < tot_uniq_mat; ++j)
	    {
		    first = uniq_map[j][0] - 1;
		    second = uniq_map[j][1] - 1;
		    exp_Xbeta2_map[first][second] = exp_Xbeta2_map_new[first][second];
	    }
	}
	else
	{
		GetRNGstate();
		x = rbinom(1, accept_prob);
		PutRNGstate();
		if (x == 1)
		{
			beta2[which_beta2] = beta2_new;
			for (j = 0; j < tot_uniq_mat; ++j)
	    {
		    first = uniq_map[j][0] - 1;
		    second = uniq_map[j][1] - 1;
		    exp_Xbeta2_map[first][second] = exp_Xbeta2_map_new[first][second];
	    }
		}
	}
}

double update_lambda_bi(double *beta1, double *beta2, double a, double b, int x_length)
{
	double lambda, beta1_abs[x_length + 1], beta2_abs[x_length + 1];
	int i;

	for (i = 0; i<x_length + 1; ++i)
	{
		beta1_abs[i] = fabs(beta1[i]);
		beta2_abs[i] = fabs(beta2[i]);
	}

	GetRNGstate();
	lambda = rgamma((double)a + 2 * (1 + x_length), 1 / (sum_bi(beta1_abs, x_length + 1) + sum_bi(beta2_abs, x_length + 1) + b));
	PutRNGstate();
	return lambda;
}

void update_u(int N, double *u, double exp_u[N], double sigma_sq_u, int ***haplo_map, int **uniq_map, double *beta1, double *beta2, int x_length, int *y1_new, int *y2_new, int tot_uniq_mat, int *num_haplo_id, double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1])
{
  double g_old, g_new, f_old,f_new,u_old, u_new, accept_prob, term1[tot_uniq_mat], term1_new[tot_uniq_mat], term2[tot_uniq_mat], term2_new[tot_uniq_mat],inprod, inprod_new,theta1, theta1_new, theta2, theta2_new, dat, dat_new, sd;
	int i, j, r, first, second, update;

 GetRNGstate();
 #pragma omp parallel for ordered private(g_old, g_new, u_old, u_new, j, r, first, second, term1, term1_new, term2, term2_new, inprod, inprod_new, theta1, theta1_new, theta2, theta2_new, dat, dat_new, sd, f_old, f_new, accept_prob)
  for (i = 0; i<N; ++i)
  {
    g_old = 0;
    g_new = 0;
  	u_old = u[i];
   #pragma omp critical(u_new)
    u_new = rnorm(u[i], sqrt(fabs(u[i])));
	  g_old = -pow(u_old, 2) / (2 * sigma_sq_u);
  	g_new = -pow(u_new, 2) / (2 * sigma_sq_u);

  	for (j = 0; j < tot_uniq_mat; ++j)
  	{
	  	first = uniq_map[j][0] - 1;
	  	second = uniq_map[j][1] - 1;
	  	term1[j] = a1_map[first][second] * exp_Xbeta1_map[first][second] * exp_u[i];
	  	term1_new[j] = a1_map[first][second] * exp_Xbeta1_map[first][second] * exp(u_new);
	  	term2[j] = a2_map[first][second] * exp_Xbeta2_map[first][second] * exp_u[i];
	  	term2_new[j] = a2_map[first][second] * exp_Xbeta2_map[first][second] * exp(u_new);
	  }

  	g_old += -0.5*(log(1 + sum_bi(term1, tot_uniq_mat)) + log(1 + sum_bi(term2, tot_uniq_mat)));
  	g_new += -0.5*(log(1 + sum_bi(term1_new, tot_uniq_mat)) + log(1 + sum_bi(term2_new, tot_uniq_mat)));

		/*Y1=0 & Y2=0, i from 1 to n_00*/
		if ((y1_new[i] == 0) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta1_new = exp_Xbeta1_map[first][second] * exp(u_new);
				theta2_new = exp_Xbeta2_map[first][second] * exp(u_new);
				dat = a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2_new + theta1_new*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=1 & Y2=0, i from n_00+1 to N2(N2=n_00+n_10)*/
		if ((y1_new[i] == 1) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta1_new = exp_Xbeta1_map[first][second] * exp(u_new);
				theta2_new = exp_Xbeta2_map[first][second] * exp(u_new);
				dat = pow(theta1, 2)* a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta1_new, 2)* a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2_new + theta1_new*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=0 & Y2=1, i from N2+1 to N3(N3=n_00+n_10+n_01)*/
		if ((y1_new[i] == 0) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta1_new = exp_Xbeta1_map[first][second] * exp(u_new);
				theta2_new = exp_Xbeta2_map[first][second] * exp(u_new);
				dat = pow(theta2, 2)* a1_map[first][second] * a2_map[first][second] / (1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow(theta2_new, 2)* a1_map[first][second] * a2_map[first][second] / (1 + theta1_new + theta2_new + theta1_new*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}

		/*Y1=1 & Y2=1, i from N3+1 to N*/
		if ((y1_new[i] == 1) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				theta1_new = exp_Xbeta1_map[first][second] * exp(u_new);
				theta2_new = exp_Xbeta2_map[first][second] * exp(u_new);
				dat = pow((theta1*theta2), 2)* a1_map[first][second] * a2_map[first][second] /(1 + theta1 + theta2 + theta1*theta2);
				dat_new = pow((theta1_new*theta2_new), 2)* a1_map[first][second] * a2_map[first][second] /(1 + theta1_new + theta2_new + theta1_new*theta2_new);
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;

			}
			g_old += log(inprod);
			g_new += log(inprod_new);

		}/*finish g_old and g_new*/

	/* f_old, f_new */
	sd=sqrt(fabs(u_old));
	f_new=-log(sd)-pow((u_old-u_new),2)/(2 * pow(sd,2));
	sd=sqrt(fabs(u_new));
	f_old=-log(sd)-pow((u_old-u_new),2)/(2 * pow(sd,2));

	accept_prob = exp(g_new - g_old + f_old -f_new);

   if (accept_prob >= 1) update = 1;
   else update = rbinom(1, accept_prob);
   if (update == 1)
   {

      u[i] = u_new;
      exp_u[i] = exp(u_new);
    }
 }
 PutRNGstate();
}

double update_sigma_sq_u(double sigma_sq_u, double *u, int N)
{
	double sigma2u, u_sq_sum = 0;
	int i;
	for (i = 0; i<N; ++i)
		u_sq_sum += u[i] * u[i];
	GetRNGstate();
	sigma2u = 1 / rgamma(((double)N - 1) / 2, 2 / (u_sq_sum));
	PutRNGstate();

	return sigma2u;

}

double update_D_bi(double *freq1, double *freq2, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, int *y1_new, int *y2_new, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map)
{
	int i, j, r, first, second, update = 0;
	double prop_D, accept_prob, g_old = 0, g_new = 0, min_f1, min_f2, min_f, delta = 0.05, lower, upper, f_old = 0, f_new = 0;
	double inprod, inprod_new, theta1, theta2, deno, dat, dat_new, term1[tot_uniq_mat], term2[tot_uniq_mat], term1_new[tot_uniq_mat], term2_new[tot_uniq_mat];
 double a1_map_new[x_length + 1][x_length + 1], a2_map_new[x_length + 1][x_length + 1];
	GetRNGstate();

	min_f1 = find_min_bi(freq1, x_length + 1);
	min_f2 = find_min_bi(freq2, x_length + 1);
	if (min_f1 < min_f2)min_f = min_f1;
	else min_f = min_f2;

	lower = D - delta;
	upper = D + delta;

	if (lower < -min_f / (1 - min_f)) lower = -min_f / (1 - min_f);
	if (upper > 1) upper = 1;

	prop_D = runif(lower, upper);

	/* calculate a1, a2, a1_new, a2_new */
	for (j = 0; j<tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		a1_map_new[first][second] = calc_a_bi(freq1, uniq_map[j], prop_D);
		a2_map_new[first][second] = calc_a_bi(freq2, uniq_map[j], prop_D);
	}

	/* g_old and g_new */
  #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term1, term1_new, term2, term2_new, theta1,theta2, deno, dat, dat_new)
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
      theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
			term1[j] = theta1 * a1_map[first][second];
			term1_new[j] = theta1 * a1_map_new[first][second];
			term2[j] = theta2 * a2_map[first][second];
			term2_new[j] = theta2 * a2_map_new[first][second];

		}
		g_old += -0.5*(log(1 + sum_bi(term1, tot_uniq_mat)) + log(1 + sum_bi(term2, tot_uniq_mat)));
		g_new += -0.5*(log(1 + sum_bi(term1_new, tot_uniq_mat)) + log(1 + sum_bi(term2_new, tot_uniq_mat)));

		if ((y1_new[i] == 0) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				deno = 1 + theta1 + theta2 + theta1*theta2;
        theta_over_deno_map[i][r] = 1 / (deno);      /* update theta_over_deno_map */
				dat = (a1_map[first][second] * a2_map[first][second]) / deno;
				dat_new = (a1_map_new[first][second] * a2_map_new[first][second]) / deno;
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);

				inprod += dat;
				inprod_new += dat_new;

			}

			g_old += log(inprod);
			g_new += log(inprod_new);

		}

		if ((y1_new[i] == 1) && (y2_new[i] == 0))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				deno = 1 + theta1 + theta2 + theta1*theta2;
        theta_over_deno_map[i][r] = pow(theta1, 2) / deno;  /* update theta_over_deno_map */
				dat = pow(theta1, 2)* (a1_map[first][second] * a2_map[first][second]) / deno;
				dat_new = pow(theta1, 2)*(a1_map_new[first][second] * a2_map_new[first][second]) / deno;

				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}


		if ((y1_new[i] == 0) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				deno = 1 + theta1 + theta2 + theta1*theta2;
        theta_over_deno_map[i][r] = pow(theta2, 2) / deno;  /* update theta_over_deno_map */

				dat = pow(theta2, 2)* (a1_map[first][second] * a2_map[first][second]) / deno;
				dat_new = pow(theta2, 2)*(a1_map_new[first][second] * a2_map_new[first][second]) / deno;
				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}


		if ((y1_new[i] == 1) && (y2_new[i] == 1))
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{

				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				theta1 = exp_Xbeta1_map[first][second] * exp_u[i];
				theta2 = exp_Xbeta2_map[first][second] * exp_u[i];
				deno = 1 + theta1 + theta2 + theta1*theta2;
        theta_over_deno_map[i][r] = pow((theta1*theta2), 2) / deno;  /* update theta_over_deno_map */

				dat = pow((theta1*theta2), 2)* (a1_map[first][second] * a2_map[first][second]) / deno;
				dat_new = pow((theta1*theta2), 2)*(a1_map_new[first][second] * a2_map_new[first][second]) / deno;

				dat = pow(dat, 0.5);
				dat_new = pow(dat_new, 0.5);
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}
	}/*finish g_old and g_new*/


	/* f_old and f_new */
	f_new = 1 / (upper - lower);

	lower = prop_D - delta;
	upper = prop_D + delta;
	if (lower < -min_f / (1 - min_f)) lower = -min_f / (1 - min_f);
	if (upper > 1) upper = 1;
	f_old = 1 / (upper - lower);
	accept_prob = exp(g_new - g_old)*f_old / f_new;

	assert(-min_f/(1-min_f)  < D);
	assert(-min_f/(1-min_f)  < prop_D);
	if (accept_prob > 1) update = 1;
	else update = rbinom(1, accept_prob);

	if (update == 1)
	{
    for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a1_map[first][second] = a1_map_new[first][second];
      a2_map[first][second] = a2_map_new[first][second];
		  }
    return prop_D;
	}
	else
		return D;
	PutRNGstate();
}


void update_freq_00(double *freq_00, double *freq_10, double *freq_01, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map)
{
	int C = 1000, update = 0;
	double prop_freq_with_last[x_length + 1], g_old, g_new, accept_prob = 0, f_old, f_new, b_old[x_length + 1], b_new[x_length + 1], term1[tot_uniq_mat], term2[tot_uniq_mat], term1_new[tot_uniq_mat], term2_new[tot_uniq_mat];
	double freq_1_new[x_length + 1], freq_2_new[x_length + 1], min_f1_old, min_f1_new, min_f2_old, min_f2_new, min_f_old, min_f_new;
  double a1_map_new[x_length + 1][x_length + 1], a2_map_new[x_length + 1][x_length + 1];
	double inprod, inprod_new, dat, dat_new;
	int i, r, j, first, second;
	GetRNGstate();

	for (i = 0; i<x_length + 1; ++i)
	{
		b_old[i] = freq_00[i] * C;
	}

	dirichlet_bi(b_old, x_length + 1, prop_freq_with_last);

	/* freq_1 and freq_2 */
	for (i = 0; i < x_length + 1; ++i)
	{
		freq_1_new[i] = (prop_freq_with_last[i] * (n00)+freq_01[i] * (n01)) / (n00 + n01);
		freq_2_new[i] = (prop_freq_with_last[i] * (n00)+freq_10[i] * (n10)) / (n00 + n10);
	}

	min_f1_old = find_min_bi(freq_1, x_length + 1);
	min_f2_old = find_min_bi(freq_2, x_length + 1);

	/* check if the constraint -min f^l_k/(1-min f^l_k) < d is satisfied */
	min_f1_new = find_min_bi(freq_1_new, x_length + 1);
	min_f2_new = find_min_bi(freq_2_new, x_length + 1);

	if (min_f1_old < min_f2_old) min_f_old = min_f1_old;
	else min_f_old = min_f2_old;

	if (min_f1_new < min_f2_new) min_f_new = min_f1_new;
	else min_f_new = min_f2_new;

	assert(-min_f_old / (1 - min_f_old)  < D);
	if (-min_f_new / (1 - min_f_new)  < D)
	{
		/* needed in acceptance prob. computation */

		for (i = 0; i<x_length + 1; ++i)
		{
			b_new[i] = prop_freq_with_last[i] * C;
			assert(b_new[i] > 0);
		}

		/* calculate a1_new, a2_new */
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
			a1_map_new[first][second] = calc_a_bi(freq_1_new, uniq_map[j], D);
			a2_map_new[first][second] = calc_a_bi(freq_2_new, uniq_map[j], D);
		}

		/* g_old and g_new */
		g_old = log(1 - min_f_old);
		g_new = log(1 - min_f_new);

    #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term1, term1_new, term2, term2_new, dat, dat_new)
		for (i = 0; i < N; ++i)
		{
			for (j = 0; j<tot_uniq_mat; ++j)
			{
				first = uniq_map[j][0] - 1;
				second = uniq_map[j][1] - 1;
				term1[j] = exp_Xbeta1_map[first][second] * exp_u[i] * a1_map[first][second];
				term1_new[j] = exp_Xbeta1_map[first][second] * exp_u[i] * a1_map_new[first][second];
				term2[j] = exp_Xbeta2_map[first][second] * exp_u[i] * a2_map[first][second];
				term2_new[j] = exp_Xbeta2_map[first][second] * exp_u[i] * a2_map_new[first][second];
			}
			g_old += -0.5*(log(1 + sum_bi(term1, tot_uniq_mat)) + log(1 + sum_bi(term2, tot_uniq_mat)));
			g_new += -0.5*(log(1 + sum_bi(term1_new, tot_uniq_mat)) + log(1 + sum_bi(term2_new, tot_uniq_mat)));

				inprod = 0;
				inprod_new = 0;
				for (r = 0; r < num_haplo_id[i]; ++r)
				{
					first = haplo_map[i][r][0] - 1;
					second = haplo_map[i][r][1] - 1;
					dat = (a1_map[first][second] * a2_map[first][second])*theta_over_deno_map[i][r];
					dat_new = (a1_map_new[first][second] * a2_map_new[first][second])*theta_over_deno_map[i][r];
					dat = pow(dat, 0.5);
					dat_new = pow(dat_new, 0.5);
					inprod += dat;
					inprod_new += dat_new;
				}
				g_old += log(inprod);
				g_new += log(inprod_new);
		} /*finish g_old and g_new*/

  /* calculate f(f_00*|f_00^(t)) = f_new and f(f_00^(t)|f_00*) = f_old */
		f_old = lgammafn(C);
		f_new = f_old;
		for (i = 0; i<x_length + 1; ++i)
		{
			f_old += -lgammafn(b_new[i]);
			f_new += -lgammafn(b_old[i]);
			f_old += (b_new[i] - 1)*log(freq_00[i]);
			f_new += (b_old[i] - 1)*log(prop_freq_with_last[i]);
		}
		accept_prob = exp(g_new - g_old + f_old - f_new);

		if (accept_prob > 1) update = 1;
		else update = rbinom(1, accept_prob);
		if (update == 1)
		{
			for (i = 0; i<x_length + 1; ++i)
      {
        freq_00[i] = prop_freq_with_last[i];
        freq_1[i] = freq_1_new[i];
        freq_2[i] = freq_2_new[i];
      }
			for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a1_map[first][second] = a1_map_new[first][second];
      a2_map[first][second] = a2_map_new[first][second];
		  }
		}
	}
	PutRNGstate();
}


void update_freq_10(double *freq_00, double *freq_10, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta2_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map)
{
	int C = 1000, update = 0;
	double prop_freq_with_last[x_length + 1], g_old, g_new, accept_prob = 0, f_old, f_new, b_old[x_length + 1], b_new[x_length + 1], term[tot_uniq_mat], term_new[tot_uniq_mat];
	double freq_2_new[x_length + 1], min_f1_old, min_f1_new, min_f2_old, min_f2_new, min_f_old, min_f_new;
	double inprod, inprod_new, dat, dat_new;
	double a2_map_new[x_length + 1][x_length + 1];
	int i, r, j, first, second;
	GetRNGstate();
	for (i = 0; i<x_length + 1; ++i)
	{
		b_old[i] = freq_10[i] * C;
	}
	dirichlet_bi(b_old, x_length + 1, prop_freq_with_last);
	for (i = 0; i < x_length + 1; ++i)
	{
		freq_2_new[i] = (freq_00[i] * (n00)+prop_freq_with_last[i] * (n10)) / (n00 + n10);
	}

	min_f1_old = find_min_bi(freq_1, x_length + 1);
	min_f2_old = find_min_bi(freq_2, x_length + 1);

	min_f1_new = min_f1_old;
	min_f2_new = find_min_bi(freq_2_new, x_length + 1);

	if (min_f1_old < min_f2_old) min_f_old = min_f1_old;
	else min_f_old = min_f2_old;

	if (min_f1_new < min_f2_new) min_f_new = min_f1_new;
	else min_f_new = min_f2_new;

	assert(-min_f_old/(1-min_f_old)  < D);
	if (-min_f_new / (1 - min_f_new)  < D)
	{
		/* needed in acceptance prob. computation */

		for (i = 0; i<x_length + 1; ++i)
		{
			b_new[i] = prop_freq_with_last[i] * C;
			assert(b_new[i] > 0);
		}
		/* calculate a1_new, a2_new */
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
			a2_map_new[first][second] = calc_a_bi(freq_2_new, uniq_map[j], D);
		}
		/* g_old and g_new */
		g_old = log(1 - min_f_old);
		g_new = log(1 - min_f_new);

    #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term, term_new, dat, dat_new)
		for (i = 0; i < N; ++i)
		{
			for (j = 0; j<tot_uniq_mat; ++j)
			{
				first = uniq_map[j][0] - 1;
				second = uniq_map[j][1] - 1;
				term[j] = exp_Xbeta2_map[first][second] * exp_u[i] * a2_map[first][second];
				term_new[j] = exp_Xbeta2_map[first][second] * exp_u[i] * a2_map_new[first][second];
			}
			g_old += -0.5*log(1 + sum_bi(term, tot_uniq_mat));
			g_new += -0.5*log(1 + sum_bi(term_new, tot_uniq_mat));

				inprod = 0;
				inprod_new = 0;
				for (r = 0; r < num_haplo_id[i]; ++r)
				{
					first = haplo_map[i][r][0] - 1;
					second = haplo_map[i][r][1] - 1;
          dat = (a1_map[first][second] * a2_map[first][second])*theta_over_deno_map[i][r];
					dat_new = (a1_map[first][second] * a2_map_new[first][second])*theta_over_deno_map[i][r];
					dat = pow(dat, 0.5);
					dat_new = pow(dat_new, 0.5);
					inprod += dat;
					inprod_new += dat_new;
				}
				g_old += log(inprod);
				g_new += log(inprod_new);
		} /*finish g_old and g_new*/

		  /* calculate f(f_00*|f_00^(t)) = f_new and f(f_00^(t)|f_00*) = f_old */
		f_old = lgammafn(C);
		f_new = f_old;
		for (i = 0; i<x_length + 1; ++i)
		{
			f_old += -lgammafn(b_new[i]);
			f_new += -lgammafn(b_old[i]);
			f_old += (b_new[i] - 1)*log(freq_10[i]);
			f_new += (b_old[i] - 1)*log(prop_freq_with_last[i]);
		}
	  accept_prob = exp(g_new - g_old + f_old - f_new);

		if (accept_prob > 1) update = 1;
		else update = rbinom(1, accept_prob);
		if (update == 1)
		{
			for (i = 0; i<x_length + 1; ++i)
      {
        freq_10[i] = prop_freq_with_last[i];
        freq_2[i] = freq_2_new[i];
      }
			for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a2_map[first][second] = a2_map_new[first][second];
		  }
		}
	}
	PutRNGstate();
}

void update_freq_01(double *freq_00, double *freq_01, double *freq_1, double *freq_2, double D, int x_length, int N, int n00, int n10, int n01, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_Xbeta1_map[x_length + 1][x_length + 1], double exp_u[N], double a1_map[x_length + 1][x_length + 1], double a2_map[x_length + 1][x_length + 1], double **theta_over_deno_map)
{
	int C = 1000, update = 0;
	double prop_freq_with_last[x_length + 1], g_old, g_new, accept_prob = 0, f_old, f_new, b_old[x_length + 1], b_new[x_length + 1], term[tot_uniq_mat], term_new[tot_uniq_mat];
	double freq_1_new[x_length + 1], min_f1_old, min_f1_new, min_f2_old, min_f2_new, min_f_old, min_f_new;
	double inprod, inprod_new, dat, dat_new;
	double a1_map_new[x_length + 1][x_length + 1];
	int i, r, j, first, second;
	GetRNGstate();

	for (i = 0; i<x_length + 1; ++i)
	{
		b_old[i] = freq_01[i] * C;
	}
	dirichlet_bi(b_old, x_length + 1, prop_freq_with_last);

	/* freq_1 and freq_2 */
	for (i = 0; i < x_length + 1; ++i)
	{
		freq_1_new[i] = (freq_00[i] * (n00)+prop_freq_with_last[i] * (n01)) / (n00 + n01);
	}

	min_f1_old = find_min_bi(freq_1, x_length + 1);
	min_f2_old = find_min_bi(freq_2, x_length + 1);

	min_f1_new = find_min_bi(freq_1_new, x_length + 1);
	min_f2_new = min_f2_old;

	/* check if the constraint -min f^l_k/(1-min f^l_k) < d is satisfied */
	if (min_f1_old < min_f2_old) min_f_old = min_f1_old;
	else min_f_old = min_f2_old;

	if (min_f1_new < min_f2_new) min_f_new = min_f1_new;
	else min_f_new = min_f2_new;

	assert(-min_f_old/(1-min_f_old)  < D);
	if (-min_f_new / (1 - min_f_new)  < D)
	{
		/* needed in acceptance prob. computation */
		for (i = 0; i<x_length + 1; ++i)
		{
			b_new[i] = prop_freq_with_last[i] * C;
			assert(b_new[i] > 0);
		}
		/* calculate a1, a2, a1_new, a2_new */
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
			a1_map_new[first][second] = calc_a_bi(freq_1_new, uniq_map[j], D);
		}
		/* g_old and g_new */
		g_old = log(1 - min_f_old);
		g_new = log(1 - min_f_new);
		#pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term, term_new, dat, dat_new)
		for (i = 0; i < N; ++i)
		{
			for (j = 0; j<tot_uniq_mat; ++j)
			{
				first = uniq_map[j][0] - 1;
				second = uniq_map[j][1] - 1;
				term[j] = exp_Xbeta1_map[first][second] * exp_u[i] * a1_map[first][second];
				term_new[j] = exp_Xbeta1_map[first][second] * exp_u[i] * a1_map_new[first][second];
			}
			g_old += -0.5*log(1 + sum_bi(term, tot_uniq_mat));
			g_new += -0.5*log(1 + sum_bi(term_new, tot_uniq_mat));

				inprod = 0;
				inprod_new = 0;
				for (r = 0; r < num_haplo_id[i]; ++r)
				{
					first = haplo_map[i][r][0] - 1;
					second = haplo_map[i][r][1] - 1;
					dat = (a1_map[first][second] * a2_map[first][second])*theta_over_deno_map[i][r];
					dat_new = (a1_map_new[first][second] * a2_map[first][second])*theta_over_deno_map[i][r];
					dat = pow(dat, 0.5);
					dat_new = pow(dat_new, 0.5);
					inprod += dat;
					inprod_new += dat_new;
				}
				g_old += log(inprod);
				g_new += log(inprod_new);
		} /*finish g_old and g_new*/

		  /* calculate f(f_00*|f_00^(t)) = f_new and f(f_00^(t)|f_00*) = f_old */
		f_old = lgammafn(C);
		f_new = f_old;

		for (i = 0; i<x_length + 1; ++i)
		{
			f_old += -lgammafn(b_new[i]);
			f_new += -lgammafn(b_old[i]);
			f_old += (b_new[i] - 1)*log(freq_01[i]);
			f_new += (b_old[i] - 1)*log(prop_freq_with_last[i]);
		}
		accept_prob = exp(g_new - g_old + f_old - f_new);

		if (accept_prob > 1) update = 1;
		else update = rbinom(1, accept_prob);
		if (update == 1)
		{
   	for (i = 0; i<x_length + 1; ++i)
      {
        freq_01[i] = prop_freq_with_last[i];
        freq_1[i] = freq_1_new[i];
      }
			for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a1_map[first][second] = a1_map_new[first][second];
		  }
	}
}

	PutRNGstate();
}

/* function to calculate a_Z in likelihood */
double calc_a_bi(double *freq, int *per_freq, double D)
{
	int i, j, k;
	double a;
	i = per_freq[0];
	j = per_freq[1];
	if (i == j){
		k = 1;
	}
	else{
		k = 0;
	}
	a = k*D*freq[i - 1] + (2 - k)*(1 - D)*freq[i - 1] * freq[j - 1];
	return a;
}

/* function to find sum_without_cov of real numbers */
double sum_bi(double *x, int n)
{
	double sum = 0.0;
	int i;

	for (i = 0; i<n; ++i)
		sum = sum + x[i];

	return sum;

}


/* function to calculate min. of an array of numbers of length n */

double find_min_bi(double *arr, int n)
{
	int i;
	double min = arr[0];
	for (i = 1; i<n; ++i)
	{
		if (min > arr[i])
			min = arr[i];
	}
	return min;
}


/* function to generate from double exponential distribution */

double gen_double_exp_bi(double mean, double SD)
{
	double x, gen_exp;

	GetRNGstate();
	x = runif(0, 1);
	gen_exp = rexp(SD / sqrt(2));
	PutRNGstate();

	if (x > 0.5)
		return gen_exp + mean;
	else
		return -gen_exp + mean;
}

/* function to generate from Dirichet distribution */

void dirichlet_bi(double *param, int dim, double *gen_sample)
{
	int i;
	double gen_gamma[dim], sum_gamma;

	GetRNGstate();
	for (i = 0; i<dim; ++i)
	{
		assert(param[i] > 0);
		gen_gamma[i] = rgamma(param[i], 1);
		if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
	}
	sum_gamma = sum_bi(gen_gamma, dim);

	for (i = 0; i<dim; ++i)
	{
		gen_sample[i] = gen_gamma[i] / sum_gamma;
	}

	PutRNGstate();
}
