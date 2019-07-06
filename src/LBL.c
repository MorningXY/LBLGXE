#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>

double calc_a_without_cov(double *freq, int *per_xz, int x_length, double D);

void update_z_without_cov(double *beta, double *freq, double D, int per_y_new_without_cov, int *per_xz, int **per_xhaplo, int num_haplo_per, int x_length);

double update_lambda_without_cov(double *beta, double a, double b, int x_length);

double calc_den_post_without_cov(double *beta, double *freq, int x_length, double D, int n_cases, int tot_hap);

void update_beta_without_cov(int **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int *y_new_without_cov, int n_cases, int tot_uniq_x, int N);

double gen_double_exp_without_cov(double mean, double SD);

void update_freq_without_cov(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz);

double update_D_without_cov(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz);

void dirichlet_without_cov(double *param, int dim, double *gen_sample);

double sum_without_cov(double *x, int n);

double find_min_without_cov(double *arr, int n);

int **uniq_x_mat_without_cov, *y_new_without_cov, total_freq_update_without_cov=0;

void mcmc_without_cov(int *x, int *tot_hap, int *y, int *N, int *num_haplo_id, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x, int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out)
{
  int i, j, k, l, m, n, ***xhaplo, **xz, x_mat[*tot_hap][*x_length], n_cases=0, which_beta, it=0, it1=0, it2=0, heat_it;

  y_new_without_cov = calloc(*N, sizeof(int));

  xhaplo = calloc(*N, sizeof(int*));
  for (i =0; i<*N; ++i)
    xhaplo[i] = calloc(num_haplo_id[i], sizeof(int*));

  for (i =0; i<*N; ++i)
    for (j = 0; j < num_haplo_id[i]; ++j)
      xhaplo[i][j] = calloc(*x_length, sizeof(int));

  /* separating x vector from R as x matrix */

  l=0;
  for (j=0; j<*x_length; ++j)
    {
      for (i=0; i<*tot_hap; ++i)
	{
	  x_mat[i][j] = x[l];
	  ++l;
	}
    }

  xz = calloc(*N, sizeof(int*));
  for (i=0; i<*N; ++i)
    xz[i] = calloc(*x_length, sizeof(int));

  /* separating haplotypes per person from the x matrix and assigning z (missing haplotype) for persons with only one compatible haplotype */

  l = 0;
  for (i =0; i<*N; ++i)
    {
      y_new_without_cov[i] = y[l];
      n_cases += y_new_without_cov[i];
      for (j = 0; j < num_haplo_id[i]; ++j)
	{
	  m = 0;
	  for (k = 0; k < *x_length; ++k)
	    {
	      xhaplo[i][j][k] = x_mat[l][m];
	      if (num_haplo_id[i]==1)
		xz[i][k] = x_mat[l][m];
	      else
		xz[i][k] = 0;
	      ++m;
	    }
	  ++l;
	}
    }

   /* separating unique x vector (haplotypes) from R as unique x matrix (to be used in denominator calculation) */

   uniq_x_mat_without_cov = calloc(*tot_uniq_x, sizeof(int*));
   for (i=0; i<*tot_uniq_x; ++i)
     uniq_x_mat_without_cov[i] = calloc(*x_length, sizeof(int));

  l=0;
  for (i=0; i<*tot_uniq_x; ++i)
    {
      for (j=0; j<*x_length; ++j)
	{
	 uniq_x_mat_without_cov[i][j] = uniq_x[l];
	  ++l;
	}
    }

  /********************** start MCMC here ***************************/
  heat_it = 0; /* tracks # of heating iterations */
  for (n=0; n<*NUM_IT; ++n)
    {
      heat_it = heat_it+1;
      if (heat_it == 101) heat_it = 1;

      /* assigning or updating z (missing haplotype) for persons with more than one compatible haplotype */

      for (i =0; i<*N; ++i)
	{
	  if (num_haplo_id[i]>1)
	    {
		  update_z_without_cov(beta, freq, *D, y_new_without_cov[i], xz[i], xhaplo[i], num_haplo_id[i], *x_length);
	    }
	}

      /* update beta parameters */

      for (i=0; i<*x_length; ++i)
	{
	  which_beta=i;
	  update_beta_without_cov(xz, beta, which_beta, *lambda, freq, *D, *x_length, y_new_without_cov, n_cases, *tot_uniq_x, *N);
	}

      /* update lambda */

      *lambda = update_lambda_without_cov(beta, *a, *b, *x_length);

       /* update frequencies and D */

      update_freq_without_cov(freq, beta, *D, *x_length, *N, n_cases, *tot_uniq_x, xz);

       /* update D parameter */

      *D = update_D_without_cov(freq, beta, *D, *x_length, *N, n_cases, *tot_uniq_x, xz);

      if (n >= *BURN_IN)
	{
	  for (i=0; i<*x_length; ++i)
	    {
	      beta_out[it] = beta[i];
	     ++it;
	    }

	  lambda_out[it2] = *lambda;

	  for (i=0; i<*x_length+1; ++i)
	    {
	      freq_out[it1] = freq[i];
	      ++it1;
	    }
	  D_out[it2] = *D;
	  ++it2;
	}

    }
}

void update_z_without_cov(double *beta, double *freq, double D, int per_y_new_without_cov, int *per_xz, int **per_xhaplo, int num_haplo_per, int x_length)
{
  int i, k;
  double prob[num_haplo_per], cum_prob[num_haplo_per], sum_without_cov_prob, x, a[num_haplo_per];

  for (i=0; i<num_haplo_per; ++i)
    {
      a[i] = calc_a_without_cov(freq, per_xhaplo[i], x_length, D);
      prob[i] = a[i];

      if (per_y_new_without_cov == 1)
	{
	  for (k=0; k < x_length; ++k)
	    {
	      prob[i] = prob[i]*exp(per_xhaplo[i][k]*beta[k]);
	    }
	}
    }

  sum_without_cov_prob = sum_without_cov(prob, num_haplo_per);

  for (i=0; i<num_haplo_per; ++i)
    {
      prob[i] = prob[i]/sum_without_cov_prob;
    }


  cum_prob[0]= prob[0];
  for (i=1; i<num_haplo_per; ++i)
    {
      cum_prob[i] = cum_prob[i-1]+prob[i];
     }

  GetRNGstate();
  x=runif(0,1);
  PutRNGstate();

  if (x < cum_prob[0])
    {
      for (k=0; k<x_length; ++k)
	per_xz[k] = per_xhaplo[0][k];
    }

  for (i=1; i<num_haplo_per; ++i)
    {
      if (x > cum_prob[i-1] && x < cum_prob[i])
	for (k=0; k<x_length; ++k)
	  per_xz[k] = per_xhaplo[i][k];
    }


}

void update_beta_without_cov(int **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int *y_new_without_cov, int n_cases, int tot_uniq_x, int N)
{
  double beta_new, g_old, g_new, f_old, f_new, beta_new_vec[x_length], SD, accept_prob;
  int i, x;

  beta_new = gen_double_exp_without_cov(beta[which_beta], sqrt(fabs(beta[which_beta])));
  g_old = -lambda*fabs(beta[which_beta]);
  g_new = -lambda*fabs(beta_new);
  for (i=0; i<N; ++i)
    {
      if (y_new_without_cov[i]==1)
	{
	  g_old += xz[i][which_beta]*beta[which_beta];
	  g_new += xz[i][which_beta]*beta_new;
	}
    }

  g_old = g_old-calc_den_post_without_cov(beta, freq, x_length, D, n_cases, tot_uniq_x);

  for (i=0; i<x_length; ++i)
    beta_new_vec[i] = beta[i];
  beta_new_vec[which_beta] = beta_new;

  g_new = g_new-calc_den_post_without_cov(beta_new_vec, freq, x_length, D, n_cases, tot_uniq_x);

  SD = sqrt(fabs(beta_new));
  f_old = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);

  SD = sqrt(fabs(beta[which_beta]));
  f_new = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
  accept_prob = exp(g_new-g_old)*f_old/f_new;

  if (accept_prob > 1)
    beta[which_beta] = beta_new;
  else
    {
      GetRNGstate();
      x = rbinom(1,accept_prob);
      PutRNGstate();
      if (x == 1) beta[which_beta] = beta_new;
    }
}

double update_lambda_without_cov(double *beta, double a, double b, int x_length)
{
  double lambda, beta_abs[x_length];
  int i;

  for (i=0; i<x_length; ++i)
    beta_abs[i] = fabs(beta[i]);
  GetRNGstate();
  lambda=rgamma((double) a+x_length, 1/(sum_without_cov(beta_abs, x_length)+b));
  PutRNGstate();

  return lambda;
}

void update_freq_without_cov(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz)
{
  int i, C=10000, update=0;
  double prop_freq_with_last[x_length+1], g_old, g_new, accept_prob=0, f_old, f_new, b_old[x_length+1], b_new[x_length+1], min_f_old, min_f_new;

  GetRNGstate();

  for (i=0; i<x_length+1; ++i)
    {
      b_old[i] = freq[i]*C;
    }

  dirichlet_without_cov(b_old, x_length+1, prop_freq_with_last);
  min_f_old = find_min_without_cov(freq, x_length+1);

  /* check if the constraint -min fk/(1-min fk) < d is satisfied */

  min_f_new = find_min_without_cov(prop_freq_with_last, x_length+1);

  assert(-min_f_old/(1-min_f_old)  < D);
  if (-min_f_new/(1-min_f_new)  < D)
    {
      /* needed in acceptance prob. computation */

      for (i=0; i<x_length+1; ++i)
	{
	  b_new[i] = prop_freq_with_last[i]*C;
	  assert(b_new[i] > 0);
	}

      /* calculate g(f^(t)) and g(f*) */

      g_old = -calc_den_post_without_cov(beta, freq, x_length, D, n_cases, tot_uniq_x)+log(1-min_f_old);

      g_new = -calc_den_post_without_cov(beta, prop_freq_with_last, x_length, D, n_cases, tot_uniq_x)+log(1-min_f_new);

      assert(calc_den_post_without_cov(beta, prop_freq_with_last, x_length, D, n_cases, tot_uniq_x) > -pow(10,10));

      for (i=0; i<N; ++i)
	{
	  g_old += log(calc_a_without_cov(freq, xz[i], x_length, D));
	  g_new += log(calc_a_without_cov(prop_freq_with_last, xz[i], x_length, D));
	}

      /* calculate f(f*|f^(t)) = f_new and f(f^(t)|f*) = f_old */

      f_old = lgammafn(C);
      f_new = f_old;

      for (i=0; i<x_length+1; ++i)
	{
	  f_old += (b_new[i]-1)*log(freq[i]) - lgammafn(b_new[i]);
	  f_new += (b_old[i]-1)*log(prop_freq_with_last[i]) - lgammafn(b_old[i]);
	}

      accept_prob = exp(g_new-g_old+f_old-f_new);

      if (accept_prob > 1) update = 1;
      else update = rbinom(1, accept_prob);
      if (update ==1)
	{
	  for (i=0; i<x_length+1; ++i)
	    freq[i] = prop_freq_with_last[i];
	  ++total_freq_update_without_cov;
	}
    }

  PutRNGstate();
}

double update_D_without_cov(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz)
{
  int i, update = 0;
  double prop_D, accept_prob, g_old, g_new, min_f, delta=0.05, lower, upper, f_old, f_new;

  GetRNGstate();

  min_f = find_min_without_cov(freq, x_length+1);
  lower = D-delta;
  upper = D+delta;

  if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
  if (upper > 1) upper = 1;

   prop_D = runif(lower, upper);

  g_old = -calc_den_post_without_cov(beta, freq, x_length, D, n_cases, tot_uniq_x);
  g_new = -calc_den_post_without_cov(beta, freq, x_length, prop_D, n_cases, tot_uniq_x);

  for (i=0; i<N; ++i)
    {
      g_old += log(calc_a_without_cov(freq, xz[i], x_length, D));
      g_new += log(calc_a_without_cov(freq, xz[i], x_length, prop_D));
    }

  f_new = 1/(upper-lower);

  lower = prop_D-delta;
  upper = prop_D+delta;
  if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
  if (upper > 1) upper = 1;

  f_old = 1/(upper-lower);

  accept_prob = exp(g_new-g_old)*f_old/f_new;
  assert(-min_f/(1-min_f)  < D);
  assert(-min_f/(1-min_f)  < prop_D);
  if (accept_prob > 1) update = 1;
  else update = rbinom(1, accept_prob);

  if (update == 1)
    return prop_D;
  else
    return D;

  PutRNGstate();
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, f, and d */

double calc_den_post_without_cov(double *beta, double *freq, int x_length, double D, int n_cases, int tot_hap)
{
  int i, k;
  double den=0, term[tot_hap];

  for (i=0; i<tot_hap; ++i)
    {
      term[i] = calc_a_without_cov(freq, uniq_x_mat_without_cov[i], x_length, D);
      for (k=0; k < x_length; ++k)
	    {
      	      term[i] = term[i]*exp(uniq_x_mat_without_cov[i][k]*beta[k]);
	    }
    }

  den = n_cases*log(sum_without_cov(term, tot_hap));
  return den;
}

/* Calculate a(F) that is in the denominator of the likelihood*/

double calc_a_without_cov(double *freq, int *per_xz, int x_length, double D)
{
  int i, j, num_haplo=0;
  double a;

  for (i=0; i < x_length; ++i)
    {
      if (per_xz[i]==1)
	{
	  a = 2*(1-D)*freq[i];
	  ++num_haplo;

	  for (j=i+1; j<x_length; ++j)
	    {
	      if (per_xz[j]==1)
		{
		  a = a*freq[j];
		  ++num_haplo;
		  break;
		}
	    }
	  if (num_haplo==2)
	    break;
	  else
	    {
	      a = a*freq[x_length];
	      ++num_haplo;
	      break;
	    }
	}

      if (per_xz[i]==2)
	{
	  a = D*freq[i]+(1-D)*pow(freq[i],2);
	  num_haplo = num_haplo+2;
	  break;
	}
    }
  if (num_haplo != 2)
    {
      a = D*freq[x_length] + (1-D)*pow(freq[x_length], 2);
    }

  return a;
}

/* function to find sum_without_cov of real numbers */

double sum_without_cov(double *x, int n)
{
  double sum_without_cov=0.0;
  int i;

  for (i=0; i<n ; ++i)
    sum_without_cov = sum_without_cov + x[i];

  return sum_without_cov;

}

/* function to calculate min. of an array of numbers of length n */

double find_min_without_cov(double *arr, int n)
{
  int i;
  double min=arr[0];
  for(i=1;i<n; ++i)
    {
      if(min > arr[i])
          min = arr[i];
    }
  return min;
}

/* function to generate from double exponential distribution */

double gen_double_exp_without_cov(double mean, double SD)
{
  double x, gen_exp;

  GetRNGstate();
  x = runif(0,1);
  gen_exp = rexp(SD/sqrt(2));
  PutRNGstate();

  if(x > 0.5)
    return gen_exp+mean;
  else
    return -gen_exp+mean;
}

/* function to generate from Dirichet distribution */

void dirichlet_without_cov(double *param, int dim, double *gen_sample)
{
  int i;
  double gen_gamma[dim], sum_without_cov_gamma;

  GetRNGstate();
  for (i=0; i<dim; ++i)
    {
      assert(param[i]>0);
      gen_gamma[i] = rgamma(param[i], 1);
      if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
    }
  sum_without_cov_gamma = sum_without_cov(gen_gamma, dim);

  for (i=0; i<dim; ++i)
    {
      gen_sample[i] = gen_gamma[i]/sum_without_cov_gamma;
    }

  PutRNGstate();
}

