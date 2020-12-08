#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <omp.h>
#include <assert.h>

void update_beta_b_bibc(int ***haplo_map, int **uniq_map, double *beta_b, int which_beta, int N, int x_length, double lambda_b, double exp_u[N], int *yb_new, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double **exp_c_ir);
double update_sigma_sq_c_bibc(int ***haplo_map, int N, int x_length, double a0, double b0, double sigma_sq_c, int *yb_new, double *yc_new, double *u, double exp_u[N], int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new);
void update_beta_c_bibc(int ***haplo_map, int **uniq_map, double *beta_c, int which_beta, int N, int x_length, double lambda_c, double exp_u[N], double sigma_sq_c, int *yb_new, double *yc_new, double *u, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new);
double update_lambda_bibc(double *beta, double a, double b, int x_length);
void update_u_bibc(int N, double *u, double exp_u[N], int *yb_new, double *yc_new, double lambda_c, double sigma_u, double sigma_sq_c, int ***haplo_map, int **uniq_map, int x_length, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new);
double update_sigma_u_bibc(double sigma_u, double *u, int N);
double update_D_bibc(double *freq, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, int *yb_new, double exp_u[N], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double a_map[x_length + 1][x_length + 1], double **exp_c_ir, double **theta_exp_map);
void update_freq_bibc(double *freq, int *yb_new, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_u[N], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double a_map[x_length + 1][x_length + 1], double **exp_c_ir, double **theta_exp_map);
double calc_a_bibc(double *freq, int *per_freq, double D);
double sum_bibc(double *x, int n);
void dirichlet_without_cov_bibc(double *param, int dim, double *gen_sample);
double find_min_without_cov_bibc(double *arr, int n);
double gen_double_exp_without_cov_bibc(double mean, double SD);



void mcmc_BiBC(int *haplotype_map, int *tot_hap, int *yb, double *yc, int *N, int *num_haplo_id, int *x_length, double *freq, double *D, double *beta_b, double *beta_c, double *a, double *b, double *a0, double *b0, int *unique_map, int *tot_uniq_mat, double *lambda_b, double *lambda_c, double *u, double *sigma_u, double *sigma_sq_c, int *NUM_IT, int *BURN_IN, double *beta_b_out, double *beta_c_out, double *lambda_b_out, double *lambda_c_out, double *freq_out, double *D_out, double *sigma_u_out, double *sigma_sq_c_out)
{
	/*  sum_indepmary of notations
	yb[i]: responce(case/control) of the disease with binary status for each subject (i<*N)
	yc[i]: responce(case/control) of the disease with continuous status for each subject (i<*N)

	*N: num of subjects
	*tot_hap: nrow of the design matrix
	num_haplo_id[i]: num of rows for each subject (i<*N)

	*x_length: length of design matrix (length of beta vector=*x_length+1)

	haplo_map[i][j][k]: map of haplotypes (i<*N, j<num_haplo_id[i], k<2)
	uniq_map[i][j]: map of haplotypes
	*tot_uniq_mat: num of unique rows of haplo.map
	*/
	int i, j, k, l, m, n, first, second, which_beta, it = 0, it1 = 0, it2 = 0/*,it3 = 0*/;
	int ***haplo_map, h_mat[*tot_hap][2], **uniq_map, yb_new[*N];
	double yc_new[*N], temp1, temp2, exp_Xbeta_b_map[*x_length + 1][*x_length + 1], Xbeta_c_map[*x_length + 1][*x_length + 1], exp_u[*N],   a_map[*x_length + 1][*x_length + 1];
  double **theta_exp_map, **exp_c_ir, **exp_c_ir_new;

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
		yb_new[i] = yb[l];
		yc_new[i] = yc[l];
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
  theta_exp_map = calloc(*N, sizeof(double*));
  exp_c_ir = calloc(*N, sizeof(double*));
  exp_c_ir_new = calloc(*N, sizeof(double*));
	for (i = 0; i<*N; ++i)
	{
     exp_u[i] = exp(u[i]);
     theta_exp_map[i] = calloc(num_haplo_id[i], sizeof(double));
     exp_c_ir[i] = calloc(num_haplo_id[i], sizeof(double));
     exp_c_ir_new[i] = calloc(num_haplo_id[i], sizeof(double));
	}

  /* initial exp_Xbeta_b_map and  Xbeta_c_map */
	for (j = 0; j < *tot_uniq_mat; ++j)
	{
     first = uniq_map[j][0] - 1;
     second = uniq_map[j][1] - 1;
     a_map[first][second] = calc_a_bibc(freq, uniq_map[j], *D);
     temp1 = beta_b[0];
	   temp2 = beta_c[0];
     if (first < *x_length  && second < *x_length)
     {
	    	temp1 += beta_b[first + 1] + beta_b[second + 1];
		    temp2 += beta_c[first + 1] + beta_c[second + 1];
     }
    else if (first == *x_length && second < *x_length)
    {
 			temp1 += beta_b[second + 1];
   	  temp2 += beta_c[second + 1];
    }
    else if (first < *x_length && second == *x_length)
    {
   	  temp1 += beta_b[first + 1];
   	  temp2 += beta_c[first + 1];
   	}
    exp_Xbeta_b_map[first][second] = exp(temp1);
    Xbeta_c_map[first][second] = temp2;
	}


 /* initial exp_c_ir */
 #pragma omp parallel for private(j, first, second)
 for (i = 0; i < *N; ++i)
 {
      for (j = 0; j < num_haplo_id[i]; ++j)
			{
				first = haplo_map[i][j][0] - 1;
				second = haplo_map[i][j][1] - 1;
        exp_c_ir[i][j] = exp(((yc_new[i]-u[i]*(sqrt(3)/M_PI)*sqrt(*sigma_sq_c))*Xbeta_c_map[first][second] - 0.5*pow(Xbeta_c_map[first][second], 2))/(*sigma_sq_c));
			}
 }

	/*---------------------start MCMC here------------------------------*/
	for (n = 0; n<*NUM_IT; ++n)
	{
    /* update beta_b parameters */
    for (i = 0; i<*x_length + 1; ++i)
		{
      which_beta = i;
      update_beta_b_bibc(haplo_map, uniq_map, beta_b, which_beta, *N, *x_length, *lambda_b, exp_u, yb_new, *tot_uniq_mat, num_haplo_id, a_map, exp_Xbeta_b_map, exp_c_ir);
		}

    /* update beta_c parameters */
		for (i = 0; i<*x_length + 1; ++i)
		{
			which_beta = i;
      update_beta_c_bibc(haplo_map, uniq_map, beta_c, which_beta, *N, *x_length, *lambda_c, exp_u, *sigma_sq_c, yb_new, yc_new, u, *tot_uniq_mat, num_haplo_id, a_map, exp_Xbeta_b_map, Xbeta_c_map, exp_c_ir, exp_c_ir_new);
		}

    /* update sigma_sq_c parameters */
    *sigma_sq_c = update_sigma_sq_c_bibc(haplo_map, *N, *x_length, *a0, *b0, *sigma_sq_c, yb_new, yc_new, u, exp_u, num_haplo_id, a_map, exp_Xbeta_b_map, Xbeta_c_map, exp_c_ir, exp_c_ir_new);


    /* update lambda_b and lambda_c */
		*lambda_b = update_lambda_bibc(beta_b, *a, *b, *x_length);
    *lambda_c = update_lambda_bibc(beta_c, *a, *b, *x_length);

    /* update vector u */
    update_u_bibc(*N, u, exp_u, yb_new, yc_new, *lambda_c, *sigma_u, *sigma_sq_c, haplo_map, uniq_map, *x_length, *tot_uniq_mat, num_haplo_id, a_map, exp_Xbeta_b_map, Xbeta_c_map, exp_c_ir, exp_c_ir_new);

		/*Update sigma_u */
		*sigma_u = update_sigma_u_bibc(*sigma_u, u, *N);

		/* update D */
		/* calculate theta_exp_map at the same time -- used in updatings of freq's and D */
    *D = update_D_bibc(freq, *D, *x_length, *N, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, yb_new, exp_u, exp_Xbeta_b_map, a_map, exp_c_ir, theta_exp_map);

		/* Updating freq */
     update_freq_bibc(freq, yb_new, *D, *x_length, *N, num_haplo_id, *tot_uniq_mat, haplo_map, uniq_map, exp_u, exp_Xbeta_b_map, a_map, exp_c_ir, theta_exp_map);

		/* outputs */
		if (n >= *BURN_IN)
		{
			for (i = 0; i<*x_length + 1; ++i)
			{
				beta_b_out[it] = beta_b[i];
				beta_c_out[it] = beta_c[i];
				++it;
			}

			lambda_b_out[it2] = *lambda_b;
      lambda_c_out[it2] = *lambda_c;
			for (i = 0; i<*x_length + 1; ++i)
			{
				freq_out[it1] = freq[i];
				++it1;
			}

			D_out[it2] = *D;
			sigma_u_out[it2] = *sigma_u;
      sigma_sq_c_out[it2] = *sigma_sq_c;
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
    free(theta_exp_map[i]);
    free(exp_c_ir[i]);
    free(exp_c_ir_new[i]);
  }
	for (i = 0; i < *tot_uniq_mat; ++i) free(uniq_map[i]);
	free(uniq_map);
  free(haplo_map);
  free(theta_exp_map);
  free(exp_c_ir);
  free(exp_c_ir_new);
}


void update_beta_b_bibc(int ***haplo_map, int **uniq_map, double *beta_b, int which_beta, int N, int x_length, double lambda_b, double exp_u[N], int *yb_new, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double **exp_c_ir)
{
	double beta_b_new, g_old = 0, g_new = 0, f_old, f_new, beta_b_new_vec[x_length + 1], SD, accept_prob, term[tot_uniq_mat], term_new[tot_uniq_mat];
	double inprod, inprod_new, temp_new, dat, dat_new;
  double exp_Xbeta_b_map_new[x_length + 1][x_length + 1];
	int i, j, r, update, first, second;

	beta_b_new = gen_double_exp_without_cov_bibc(beta_b[which_beta], sqrt(fabs(beta_b[which_beta])));

	for (i = 0; i<x_length + 1; ++i)
		beta_b_new_vec[i] = beta_b[i];
	beta_b_new_vec[which_beta] = beta_b_new;

	/* calculate exp_Xbeta_b_map_new */
  #pragma omp parallel for private(first, second, temp_new)
	for (j = 0; j < tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		temp_new = beta_b_new_vec[0];
		if (first < x_length  && second < x_length)
		{
			temp_new += beta_b_new_vec[first + 1] + beta_b_new_vec[second + 1];
		}
		else if (first == x_length && second < x_length)
		{
			temp_new += beta_b_new_vec[second + 1];
		}
		else if (first < x_length && second == x_length)
		{
			temp_new += beta_b_new_vec[first + 1];
		}
		exp_Xbeta_b_map_new[first][second] = exp(temp_new); /* update exp_Xbeta_b_map_new here, already have exp_Xbeta_b_map and Xbeta_c_map */
	}

	/*g_old and g_new*/
	g_old = -lambda_b*fabs(beta_b[which_beta]);
	g_new = -lambda_b*fabs(beta_b_new);

  #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term, term_new, dat, dat_new)
	for (i = 0; i < N; ++i)
	{
    for (j = 0; j < tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;

			term[j] = a_map[first][second] *  exp_Xbeta_b_map[first][second] * exp_u[i];
			term_new[j] = a_map[first][second] * exp_Xbeta_b_map_new[first][second] * exp_u[i];
		}
		g_old += -log(1 + sum_bibc(term, tot_uniq_mat));
		g_new += -log(1 + sum_bibc(term_new, tot_uniq_mat));

		if (yb_new[i] == 1)
		{
			inprod = 0;
			inprod_new = 0;
			for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
				dat = a_map[first][second] * exp_Xbeta_b_map[first][second]*exp_u[i] * exp_c_ir[i][r];
				dat_new = a_map[first][second] * exp_Xbeta_b_map_new[first][second]*exp_u[i] * exp_c_ir[i][r];
				inprod += dat;
				inprod_new += dat_new;
			}
			g_old += log(inprod);
			g_new += log(inprod_new);
		}
	} /*finish g_old and g_new*/

	  /*f_old and f_new*/
	SD = sqrt(fabs(beta_b_new));
	f_old = exp(-sqrt(2)*fabs(beta_b[which_beta] - beta_b_new) / SD) / (sqrt(2)*SD);
	SD = sqrt(fabs(beta_b[which_beta]));
	f_new = exp(-sqrt(2)*fabs(beta_b[which_beta] - beta_b_new) / SD) / (sqrt(2)*SD);
	accept_prob = exp(g_new - g_old)*(f_old / f_new);

  if (accept_prob > 1)
	{
		update = 1;
	}
	else
	{
		GetRNGstate();
		update = rbinom(1, accept_prob);
		PutRNGstate();
   }
	if (update == 1)
	{
	  beta_b[which_beta] = beta_b_new;
    for (j = 0; j < tot_uniq_mat; ++j)
    {
      first = uniq_map[j][0] - 1;
      second = uniq_map[j][1] - 1;
      exp_Xbeta_b_map[first][second] = exp_Xbeta_b_map_new[first][second];
    }
	}
}

double update_sigma_sq_c_bibc(int ***haplo_map, int N, int x_length, double a0, double b0, double sigma_sq_c, int *yb_new, double *yc_new, double *u, double exp_u[N], int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new)
{
	double sigma_sq_c_new, g_old = 0, g_new = 0, accept_prob, inprod, inprod_new, temp, dat, dat_new, sd, f_new, f_old;
	int i, r, update, first, second;

  GetRNGstate();
  sigma_sq_c_new = rnorm(sigma_sq_c, fabs(sigma_sq_c)/15);
  PutRNGstate();
  if (sigma_sq_c_new <= 0) return sigma_sq_c;

 g_old = -(N+1)*log(sigma_sq_c)/2;
 g_new = -(N+1)*log(sigma_sq_c_new)/2;

  #pragma omp parallel for reduction(+:g_old, g_new) private(r, first, second, inprod, inprod_new, temp, dat, dat_new)
	for (i = 0; i < N; ++i)
	{
    g_old += -pow((yc_new[i]-u[i]*(sqrt(3)/M_PI)*sqrt(sigma_sq_c)), 2)/(2 * sigma_sq_c);
		g_new += -pow((yc_new[i]-u[i]*(sqrt(3)/M_PI)*sqrt(sigma_sq_c_new)), 2)/(2 * sigma_sq_c_new);

    inprod = 0;
		inprod_new = 0;
     for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
        temp = exp(((yc_new[i]-u[i]*(sqrt(3)/M_PI)*sqrt(sigma_sq_c_new))*Xbeta_c_map[first][second] - 0.5*pow(Xbeta_c_map[first][second], 2))/sigma_sq_c_new);
        exp_c_ir_new[i][r] = temp;
        dat = a_map[first][second] * exp_c_ir[i][r];
		    dat_new = a_map[first][second] * temp;
				if (yb_new[i] == 1)
		    {
			    dat = dat * exp_Xbeta_b_map[first][second] * exp_u[i];
				  dat_new = dat_new * exp_Xbeta_b_map[first][second] * exp_u[i];
		    }
				inprod += dat;
				inprod_new += dat_new;
			}
    g_old += log(inprod);
		g_new += log(inprod_new);
	} /*finish g_old and g_new*/

  sd=fabs(sigma_sq_c)/15;
	f_new=-log(sd)-pow((sigma_sq_c-sigma_sq_c_new),2)/(2 * pow(sd,2));
	sd=fabs(sigma_sq_c_new)/15;
	f_old=-log(sd)-pow((sigma_sq_c-sigma_sq_c_new),2)/(2 * pow(sd,2));
	accept_prob = exp(g_new - g_old + f_old -f_new);

	if (accept_prob > 1)
	{
		update = 1;
	}
	else
	{
		GetRNGstate();
		update = rbinom(1, accept_prob);
		PutRNGstate();
  }
	if (update == 1){
    for (i = 0; i < N; ++i){
      for (r = 0; r < num_haplo_id[i]; ++r){
        exp_c_ir[i][r] = exp_c_ir_new[i][r];
      }
    }
    return sigma_sq_c_new;
	}
 else{
   return sigma_sq_c;
 }
}

void update_beta_c_bibc(int ***haplo_map, int **uniq_map, double *beta_c, int which_beta, int N, int x_length, double lambda_c, double exp_u[N], double sigma_sq_c, int *yb_new, double *yc_new, double *u, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new)
{
	double beta_c_new, g_old = 0, g_new = 0, f_old, f_new, beta_c_new_vec[x_length + 1], SD, accept_prob, inprod, inprod_new, temp_new, temp1_new, Xbeta_c_map_new[x_length + 1][x_length + 1];
	int i, j, r, update, first, second;

	beta_c_new = gen_double_exp_without_cov_bibc(beta_c[which_beta], sqrt(fabs(beta_c[which_beta])));

	for (i = 0; i<x_length + 1; ++i)
		beta_c_new_vec[i] = beta_c[i];
	beta_c_new_vec[which_beta] = beta_c_new;

	/*calculate Xbeta_b_map_new*/
  #pragma omp parallel for private(first, second, temp_new)
	for (j = 0; j < tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		temp_new = beta_c_new_vec[0];
		if (first < x_length  && second < x_length)
		{
			temp_new += beta_c_new_vec[first + 1] + beta_c_new_vec[second + 1];
		}
		else if (first == x_length && second < x_length)
		{
			temp_new += beta_c_new_vec[second + 1];
		}
		else if (first < x_length && second == x_length)
		{
			temp_new += beta_c_new_vec[first + 1];
		}
		Xbeta_c_map_new[first][second] = temp_new; /* update Xbeta_c_map_new here, already have Xbeta_c_map */
	}

	/*g_old and g_new*/
	g_old = -lambda_c*fabs(beta_c[which_beta]);
	g_new = -lambda_c*fabs(beta_c_new);

  #pragma omp parallel for reduction(+:g_old, g_new) private(r, first, second, inprod, inprod_new, temp1_new)
	for (i = 0; i < N; ++i)
	{
    inprod = 0;
		inprod_new = 0;
   for (r = 0; r < num_haplo_id[i]; ++r)
	{
	   first = haplo_map[i][r][0] - 1;
	   second = haplo_map[i][r][1] - 1;
     temp1_new = exp(((yc_new[i]-u[i]*(sqrt(3)/M_PI)*sqrt(sigma_sq_c))*Xbeta_c_map_new[first][second] - 0.5*pow(Xbeta_c_map_new[first][second], 2))/sigma_sq_c);
     exp_c_ir_new[i][r] = temp1_new;
     if (yb_new[i] == 0)
     {
			  inprod += a_map[first][second] * exp_c_ir[i][r];
			  inprod_new += a_map[first][second] * temp1_new;
     }
     if (yb_new[i] == 1)
     {
			  inprod += a_map[first][second] * exp_Xbeta_b_map[first][second] * exp_u[i] * exp_c_ir[i][r];
			  inprod_new += a_map[first][second] * exp_Xbeta_b_map[first][second] * exp_u[i] * temp1_new;
     }
	}
			g_old += log(inprod);
			g_new += log(inprod_new);
	} /*finish g_old and g_new*/

	  /*f_old and f_new*/
	SD = sqrt(fabs(beta_c_new));
	f_old = exp(-sqrt(2)*fabs(beta_c[which_beta] - beta_c_new) / SD) / (sqrt(2)*SD);
	SD = sqrt(fabs(beta_c[which_beta]));
	f_new = exp(-sqrt(2)*fabs(beta_c[which_beta] - beta_c_new) / SD) / (sqrt(2)*SD);
	accept_prob = exp(g_new - g_old)*(f_old / f_new);

	if (accept_prob > 1)
	{
		update = 1;
	}
	else
	{
		GetRNGstate();
		update = rbinom(1, accept_prob);
		PutRNGstate();
  }
	if (update == 1)
	{
		beta_c[which_beta] = beta_c_new;
    for (j = 0; j < tot_uniq_mat; ++j)
    {
      first = uniq_map[j][0] - 1;
	    second = uniq_map[j][1] - 1;
	    Xbeta_c_map[first][second] = Xbeta_c_map_new[first][second];
    }
    for (i = 0; i < N; ++i)
    {
      for (r = 0; r < num_haplo_id[i]; ++r)
      {
        exp_c_ir[i][r] = exp_c_ir_new[i][r];
      }
    }
	}
}



double update_lambda_bibc(double *beta, double a, double b, int x_length)
{
	double lambda, beta_abs[x_length + 1];
	int i;
	for (i = 0; i<x_length + 1; ++i)
	{
		beta_abs[i] = fabs(beta[i]);
	}
	GetRNGstate();
	lambda = rgamma((double)a + 1 + x_length, 1 / (sum_bibc(beta_abs, x_length + 1) + b));
	PutRNGstate();
	return lambda;
}



void update_u_bibc(int N, double *u, double exp_u[N], int *yb_new, double *yc_new, double lambda_c, double sigma_u, double sigma_sq_c, int ***haplo_map, int **uniq_map, int x_length, int tot_uniq_mat, int *num_haplo_id, double a_map[x_length + 1][x_length + 1], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double Xbeta_c_map[x_length + 1][x_length + 1], double **exp_c_ir, double **exp_c_ir_new)
{
  double g_old, g_new, f_old,f_new,u_old, u_new, accept_prob, term[tot_uniq_mat], term_new[tot_uniq_mat],temp, inprod, inprod_new, sd;
	int i, j, r, first, second, update;
 GetRNGstate();
 #pragma omp parallel for private(g_old, g_new, u_old, u_new, j, r, first, second, term, term_new, temp, inprod, inprod_new, sd, f_old, f_new, accept_prob)
  for (i = 0; i<N; ++i)
  {
    g_old=0;
    g_new=0;
  	u_old = u[i];
    #pragma omp critical(u_new)
    u_new = rnorm(u_old, sqrt(fabs(u_old)));
	  g_old = -pow(u_old, 2) / (2 * pow(sigma_u,2)) - pow((yc_new[i]-u_old*(sqrt(3)/M_PI)*sqrt(sigma_sq_c)), 2)/(2 * sigma_sq_c);
  	g_new = -pow(u_new, 2) / (2 * pow(sigma_u,2)) - pow((yc_new[i]-u_new*(sqrt(3)/M_PI)*sqrt(sigma_sq_c)), 2)/(2 * sigma_sq_c);

  	for (j = 0; j < tot_uniq_mat; ++j)
  	{
	  	first = uniq_map[j][0] - 1;
	  	second = uniq_map[j][1] - 1;
	  	term[j] = a_map[first][second] * exp_Xbeta_b_map[first][second] * exp_u[i];
	  	term_new[j] = a_map[first][second] * exp_Xbeta_b_map[first][second] * exp(u_new);
	  }
  	g_old -= log(1 + sum_bibc(term, tot_uniq_mat));
  	g_new -= log(1 + sum_bibc(term_new, tot_uniq_mat));

    inprod = 0;
    inprod_new = 0;
    for (r = 0; r < num_haplo_id[i]; ++r)
    {
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
        temp = exp(((yc_new[i]-u_new*(sqrt(3)/M_PI)*sqrt(sigma_sq_c))*Xbeta_c_map[first][second] - 0.5*pow(Xbeta_c_map[first][second], 2))/sigma_sq_c);
        exp_c_ir_new[i][r] = temp;
        if (yb_new[i] == 0)
		    {
			  	inprod += a_map[first][second] * exp_c_ir[i][r];
				  inprod_new += a_map[first][second] * temp;
        }
        if (yb_new[i] == 1)
		    {
			  	inprod += a_map[first][second] * exp_Xbeta_b_map[first][second] * exp_u[i] * exp_c_ir[i][r];
			  	inprod_new += a_map[first][second] * exp_Xbeta_b_map[first][second] * exp(u_new) * temp;
        }
    }
    g_old += log(inprod);
    g_new += log(inprod_new);
		/*finish g_old and g_new*/
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
      for (r = 0; r < num_haplo_id[i]; ++r)
			{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
        exp_c_ir[i][r]= exp_c_ir_new[i][r];
			}
    }
 }
 PutRNGstate();
}

double update_sigma_u_bibc(double sigma_u, double *u, int N)
{
  double sigma_u_new, u_sq_sum = 0, g_old, g_new, f_old, f_new, sd, accept_prob, update=0;
	int i;
  GetRNGstate();
  sigma_u_new = rnorm(sigma_u, fabs(sigma_u)/20);
  PutRNGstate();

  if (sigma_u_new <= 0) return sigma_u;

 for (i = 0; i<N; ++i){
 u_sq_sum += u[i] * u[i];
 }
  g_old = -N*log(sigma_u) - log(1 + pow(sigma_u/10, 2)) - u_sq_sum/(2*pow(sigma_u, 2));
  g_new = -N*log(sigma_u_new) - log(1 + pow(sigma_u_new/10, 2)) - u_sq_sum/(2*pow(sigma_u_new, 2));

	sd=fabs(sigma_u)/20;
	f_new=-log(sd)-pow((sigma_u-sigma_u_new),2)/(2 * pow(sd,2));
	sd=fabs(sigma_u_new)/20;
	f_old=-log(sd)-pow((sigma_u-sigma_u_new),2)/(2 * pow(sd,2));
	accept_prob = exp(g_new - g_old + f_old -f_new);

	if (accept_prob > 1)
	{
		update = 1;
	}
	else
	{
		GetRNGstate();
		update = rbinom(1, accept_prob);
		PutRNGstate();
  }
	if (update == 1){
    return sigma_u_new;
	}
 else{
   return sigma_u;
 }

}

double update_D_bibc(double *freq, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, int *yb_new, double exp_u[N], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double a_map[x_length + 1][x_length + 1], double **exp_c_ir, double **theta_exp_map)
{
	int i, j, r, first, second, update = 0;
	double prop_D, accept_prob, g_old = 0, g_new = 0, min_f, delta = 0.05, lower, upper, f_old = 0, f_new = 0;
	double inprod, inprod_new, theta, term[tot_uniq_mat], term_new[tot_uniq_mat];
 double a_map_new[x_length + 1][x_length + 1];
	GetRNGstate();
	min_f = find_min_without_cov_bibc(freq, x_length + 1);
	lower = D - delta;
	upper = D + delta;
	if (lower < -min_f / (1 - min_f)) lower = -min_f / (1 - min_f);
	if (upper > 1) upper = 1;
	prop_D = runif(lower, upper);

	for (j = 0; j<tot_uniq_mat; ++j)
	{
		first = uniq_map[j][0] - 1;
		second = uniq_map[j][1] - 1;
		a_map_new[first][second] = calc_a_bibc(freq, uniq_map[j], prop_D);
	}
	/* g_old and g_new */
  #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, theta, inprod, inprod_new, term, term_new)
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      theta = exp_Xbeta_b_map[first][second] * exp_u[i];
			term[j] = theta * a_map[first][second];
			term_new[j] = theta * a_map_new[first][second];
		}
		g_old -= log(1 + sum_bibc(term, tot_uniq_mat));
		g_new -= log(1 + sum_bibc(term_new, tot_uniq_mat));

    inprod = 0;
		inprod_new = 0;
    for (r = 0; r < num_haplo_id[i]; ++r)
		{
				first = haplo_map[i][r][0] - 1;
				second = haplo_map[i][r][1] - 1;
        if (yb_new[i] == 0)
        {
				  inprod += a_map[first][second] * exp_c_ir[i][r];
				  inprod_new += a_map_new[first][second] * exp_c_ir[i][r];
        }
        if (yb_new[i] == 1)
		    {
          theta_exp_map[i][r] = exp_Xbeta_b_map[first][second] * exp_u[i] * exp_c_ir[i][r];
				  inprod += a_map[first][second] * theta_exp_map[i][r];
				  inprod_new += a_map_new[first][second] * theta_exp_map[i][r];
        }
		}
   g_old += log(inprod);
   g_new += log(inprod_new);
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
	PutRNGstate();
	if (update == 1)
	{
    for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a_map[first][second] = a_map_new[first][second];
		  }
    return prop_D;
	}
	else
		return D;
}

void update_freq_bibc(double *freq, int *yb_new, double D, int x_length, int N, int *num_haplo_id, int tot_uniq_mat, int ***haplo_map, int **uniq_map, double exp_u[N], double exp_Xbeta_b_map[x_length + 1][x_length + 1], double a_map[x_length + 1][x_length + 1], double **exp_c_ir, double **theta_exp_map)
{
	int C = 10000, update = 0;
	double freq_new[x_length + 1], g_old, g_new, accept_prob = 0, f_old, f_new, b_old[x_length + 1], b_new[x_length + 1], term[tot_uniq_mat], term_new[tot_uniq_mat], min_f_old, min_f_new, a_map_new[x_length + 1][x_length + 1];
	double inprod, inprod_new;
	int i, r, j, first, second;
	GetRNGstate();

	for (i = 0; i<x_length + 1; ++i)
	{
		b_old[i] = freq[i] * C;
	}

	dirichlet_without_cov_bibc(b_old, x_length + 1, freq_new);
	min_f_old = find_min_without_cov_bibc(freq, x_length + 1);
	min_f_new = find_min_without_cov_bibc(freq_new, x_length + 1);

	if (-min_f_new / (1 - min_f_new)  < D)
	{
		/* needed in acceptance prob. computation */
		for (i = 0; i<x_length + 1; ++i)
		{
			b_new[i] = freq_new[i] * C;
		  assert(b_new[i] > 0);
		}
		/* calculate a_new */
		for (j = 0; j<tot_uniq_mat; ++j)
		{
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
			a_map_new[first][second] = calc_a_bibc(freq_new, uniq_map[j], D);
		}
		/* g_old and g_new */
		g_old = log(1 - min_f_old);
		g_new = log(1 - min_f_new);

    #pragma omp parallel for reduction(+:g_old, g_new) private(j, r, first, second, inprod, inprod_new, term, term_new)
		for (i = 0; i < N; ++i)
		{
			for (j = 0; j<tot_uniq_mat; ++j)
			{
				first = uniq_map[j][0] - 1;
				second = uniq_map[j][1] - 1;
				term[j] = exp_Xbeta_b_map[first][second] * exp_u[i] * a_map[first][second];
				term_new[j] = exp_Xbeta_b_map[first][second] * exp_u[i] * a_map_new[first][second];
			}
			g_old -= log(1 + sum_bibc(term, tot_uniq_mat));
			g_new -= log(1 + sum_bibc(term_new, tot_uniq_mat));

				inprod = 0;
				inprod_new = 0;
				for (r = 0; r < num_haplo_id[i]; ++r)
				{
					first = haplo_map[i][r][0] - 1;
					second = haplo_map[i][r][1] - 1;
					if (yb_new[i] == 0)
          {
            inprod += a_map[first][second] * exp_c_ir[i][r];
				    inprod_new += a_map_new[first][second] * exp_c_ir[i][r];
          }
		      if (yb_new[i] == 1)
          {
            inprod += a_map[first][second] * theta_exp_map[i][r];
				    inprod_new += a_map_new[first][second] * theta_exp_map[i][r];
          }
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
			f_old += (b_new[i] - 1)*log(freq[i]);
			f_new += (b_old[i] - 1)*log(freq_new[i]);
		}
		accept_prob = exp(g_new - g_old + f_old - f_new);

		if (accept_prob > 1) update = 1;
		else update = rbinom(1, accept_prob);
		if (update == 1)
		{
			for (i = 0; i<x_length + 1; ++i)
      {
        freq[i] = freq_new[i];
      }
			for (j = 0; j<tot_uniq_mat; ++j)
		  {
			first = uniq_map[j][0] - 1;
			second = uniq_map[j][1] - 1;
      a_map[first][second] = a_map_new[first][second];
		  }
		}
	}

	PutRNGstate();
}

/* function to calculate a_Z in likelihood */
double calc_a_bibc(double *freq, int *per_freq, double D)
{
	int i, j, k;
	double a;
	i = per_freq[0];
	j = per_freq[1];
	if (i == j)
	{
		k = 1;
	}
	else
	{
		k = 0;
	}
	a = k*D*freq[i - 1] + (2 - k)*(1 - D)*freq[i - 1] * freq[j - 1];
	return a;
}

/* function to find sum_without_cov of real numbers */
double sum_bibc(double *x, int n)
{
	double sum = 0.0;
	int i;

	for (i = 0; i<n; ++i)
		sum = sum + x[i];

	return sum;

}


/* function to calculate min. of an array of numbers of length n */
double find_min_without_cov_bibc(double *arr, int n)
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

double gen_double_exp_without_cov_bibc(double mean, double SD)
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

void dirichlet_without_cov_bibc(double *param, int dim, double *gen_sample)
{
	int i;
	double gen_gamma[dim], sum_without_cov_gamma;

	GetRNGstate();
	for (i = 0; i<dim; ++i)
	{
	  assert(param[i] > 0);
	  gen_gamma[i] = rgamma(param[i], 1);
	  if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
	}
	sum_without_cov_gamma = sum_bibc(gen_gamma, dim);

	for (i = 0; i<dim; ++i)
	{
		gen_sample[i] = gen_gamma[i] / sum_without_cov_gamma;
	}

	PutRNGstate();
}


