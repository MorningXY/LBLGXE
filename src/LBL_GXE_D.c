#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>

double calc_a_dep(double *freq, int *per_freq, double D);

double calc_den_post_dep(double *beta, double **freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, double a[h_length][h_length][len_E], double theta[h_length][h_length][len_E]);

double update_lambda_dep(double *beta, double a, double b, int x_length);

void update_beta_dep(double *beta, int which_beta, double lambda, double **freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int *beta_accept_count, int count);

void update_b_par_dep(double **b_par, int which_b_first, int which_b_second, double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int **b_par_accept_count, int count);

double update_D_dep(double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int len_dep, int count);

double gen_double_exp_dep(double mean, double SD);

void dirichlet_dep(double *param, int dim, double *gen_sample);

double sum_dep(double *x, int n);

double find_min_dep(double *arr, int n);

void mcmc_dep(int *y, int *N, int *n_tot, int *num_haplo_id, int *x_length, int *h_length, int *haplotype_map, int *tot_uniq_mat, int *unique_map, int *index_E, int *value_Environment, int *num_E, int *len_E, int *len_dummy, int *index_Etodep, int *value_dependent, int *len_dep, int *len_dummy_dep, double *beta, double *lambda, double *freq_vec, double *D, double *a, double *b, double *beta_out, double *lambda_out, double *b_out, double *freq_out, double *D_out, int *NUM_IT, int *BURN_IN, double *start_gamma)
{

/*  sum_depmary of notations
    y[i]: responce(case/control) for each subject (i<*N)
    *N: num of subjects
    *n_tot: nrow of the design matrix
    num_haplo_id[i]: num of rows for each subject (i<*N)

    *x_length: length of design matrix (length of beta vector=*x_length+1)
    *h_length: num of haplotypes

    haplo_map[i][j][k]: map of haplotypes (i<*N, j<num_haplo_id[i], k<2)
    uniq_map[i][j]: map of haplotypes
    *tot_uniq_mat: num of unique rows of haplo.mat

    index_E[i]: index of ith subject's E in value_E (i<*N)
    value_E[i][j]: all unique rows of cov.mat (i<*len_E, j<*len_dummy)
    *len_E: num of unique rows of cov.mat

    index_Etodep[i]: transform index_E to "index_dep" (i<*len_E)
    value_dep[i][j]: all unique rows of dep_mat (i<*len_dep, j<*len_dummy_dep)
    *len_dep: num of unique rows of dep_mat
*/

    int i, j, k, l, m, n, ***haplo_map, h_mat[*n_tot][2], **uniq_map, **value_E, **value_dep, n_cases=0, which_beta, which_b_first, which_b_second, it1=0, it2=0, it3=0, it4=0, heat_it;

    double **b_par, **freq, *den_freq;

    int count=0, *beta_accept_count, **b_par_accept_count;

    for (i =0; i<*N; ++i) n_cases += y[i];

    /* separating haplotype_map vector from R as h_mat matrix */

    l=0;
    for (j=0; j<2; ++j)
    {
        for (i=0; i<*n_tot; ++i)
        {
            h_mat[i][j] = haplotype_map[l];
            ++l;
        }
    }

    /* separating h_mat as haplo_map */

    haplo_map = calloc(*N, sizeof(int*));
    for (i =0; i<*N; ++i)
    {
        haplo_map[i] = calloc(num_haplo_id[i], sizeof(int*));
    }

    for (i =0; i<*N; ++i)
    {
        for (j = 0; j < num_haplo_id[i]; ++j)
        {
            haplo_map[i][j] = calloc(2, sizeof(int));
        }
    }

    l = 0;
    for (i =0; i<*N; ++i)
    {
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
    for (i=0; i<*tot_uniq_mat; ++i)
    {
        uniq_map[i] = calloc(2, sizeof(int));
    }

    l=0;
    for (i=0; i<*tot_uniq_mat; ++i)
    {
        for (j=0; j<2; ++j)
        {
            uniq_map[i][j] = unique_map[l];
            ++l;
        }
    }

    /* separating value_Environment vector from R as value_E matrix */

    value_E = calloc(*len_E, sizeof(int*));
    for (i=0; i<*len_E; ++i)
    {
        value_E[i] = calloc(*len_dummy, sizeof(int));
    }

    l=0;
    for (i=0; i<*len_E; ++i)
    {
        for (j=0; j<*len_dummy; ++j)
        {
            value_E[i][j] = value_Environment[l];
            ++l;
        }
    }

    /* separating value_dependent vector from R as value_dep matrix */

    value_dep = calloc(*len_dep, sizeof(int*));
    for (i=0; i<*len_dep; ++i)
    {
        value_dep[i] = calloc(*len_dummy_dep, sizeof(int));
    }

    l=0;
    for (i=0; i<*len_dep; ++i)
    {
        for (j=0; j<*len_dummy_dep; ++j)
        {
            value_dep[i][j] = value_dependent[l];
            ++l;
        }
    }

    /* get the starting point of b_par */

    b_par=calloc(*h_length-1, sizeof(double*));
    for(i=0;i<*h_length-1;++i)
    {
        b_par[i]=calloc(*len_dummy_dep+1, sizeof(double));
        b_par[i][0]=log(freq_vec[i]/freq_vec[*h_length-1]);
        for(j=1;j<*len_dummy_dep+1;++j)
            b_par[i][j]=*start_gamma;
    }

    /* get the corresponding freq from the above */

    freq=calloc(*len_dep, sizeof(double*));
    for(i=0;i<*len_dep;++i)
    {
        freq[i]=calloc(*h_length, sizeof(double));
        for(j=0;j<*h_length-1;++j)
        {
            freq[i][j]=exp(b_par[j][0]);
            for(k=1;k<*len_dummy_dep+1;++k)
                freq[i][j]*=exp(b_par[j][k]*value_dep[i][k-1]);
        }
        freq[i][*h_length-1]=1;
    }

    den_freq=calloc(*len_dep, sizeof(double));
    for(i=0;i<*len_dep;++i)
    {
        den_freq[i]=sum_dep(freq[i],*h_length);
        for(j=0;j<*h_length;++j)
            freq[i][j]/=den_freq[i];
    }

    /* prepare arrays for accept prob */

    b_par_accept_count=calloc(*h_length-1, sizeof(int*));
    for(i=0;i<*h_length-1;++i)
    {
        b_par_accept_count[i]=calloc(*len_dummy_dep+1, sizeof(int*));
        for(j=0;j<*len_dummy_dep+1;++j)
            b_par_accept_count[i][j]=0;
    }

    beta_accept_count=calloc(*x_length+1, sizeof(int));
    for(i=0;i<*x_length+1;++i)
    {
        beta_accept_count[i]=0;
    }

    /********************** start MCMC here ***************************/

    heat_it = 0; /* tracks # of heating iterations */

    for (n=0; n<*NUM_IT; ++n)
    {
        heat_it = heat_it+1;
        if (heat_it == 101) heat_it = 1;

        if (n >= *BURN_IN) {count = 1;}

        /* update beta parameters */

        for (i=0; i<*x_length+1; ++i)
        {
            which_beta=i;
            update_beta_dep(beta, which_beta, *lambda, freq, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, beta_accept_count, count);
        }

        /* update lambda */

        *lambda = update_lambda_dep(beta, *a, *b, *x_length);

        /* update b and frequencies */

        for (i=0; i<*h_length-1; ++i)
        {
            for (j=0; j<*len_dummy_dep+1; ++j)
            {
                which_b_first=i;
                which_b_second=j;
                update_b_par_dep(b_par, which_b_first, which_b_second, freq, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, value_dep, *len_dep, *len_dummy_dep, b_par_accept_count, count);
            }
        }

        /* update D parameter */

        *D = update_D_dep(freq, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, *len_dep, count);

        if (n >= *BURN_IN)
        {
            for (i=0; i<*x_length+1; ++i)
            {
                beta_out[it1] = beta[i];
                ++it1;
            }

            lambda_out[it2] = *lambda;

            for (i=0; i<*h_length-1; ++i)
            {
                for (j=0; j<*len_dummy_dep+1; ++j)
                {
                    b_out[it3] = b_par[i][j];
                    ++it3;
                }
            }

            for (i=0; i<*len_dep; ++i)
            {
                for (j=0; j<*h_length; ++j)
                {
                    freq_out[it4] = freq[i][j];
                    ++it4;
                }
            }

            D_out[it2] = *D;
            ++it2;
        }
    }
}

void update_beta_dep(double *beta, int which_beta, double lambda, double **freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int *beta_accept_count, int count)
{
    int i, j, update, first, second, third;
    double a[h_length][h_length][len_E], theta_old[h_length][h_length][len_E], theta_new[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double beta_new, beta_new_vec[x_length+1], den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;

    beta_new = gen_double_exp_dep(beta[which_beta], sqrt(fabs(beta[which_beta])));
    for (i=0; i<x_length+1; ++i)
        beta_new_vec[i] = beta[i];
    beta_new_vec[which_beta] = beta_new;

    temp_old=calloc(N, sizeof(double*));
    temp_new=calloc(N, sizeof(double*));
    for(i=0;i<N;++i)
    {
        temp_old[i]=calloc(num_haplo_id[i], sizeof(double));
        temp_new[i]=calloc(num_haplo_id[i], sizeof(double));
    }

    g_ratio=1;
    den_old=calc_den_post_dep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a, theta_old);
    den_new=calc_den_post_dep(beta_new_vec, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a, theta_new);
    for(i=0;i<N;++i)
    {
        t_old[i]=0;
        t_new[i]=0;
        for(j=0;j<num_haplo_id[i];++j)
        {
            first=haplo_map[i][j][0]-1;
            second=haplo_map[i][j][1]-1;
            third=index_E[i]-1;
            temp_old[i][j]=a[first][second][third];
            temp_new[i][j]=temp_old[i][j];
            if (y[i]==1)
            {
                temp_old[i][j]*=theta_old[first][second][third];
                temp_new[i][j]*=theta_new[first][second][third];
            }
            t_old[i]+=temp_old[i][j];
            t_new[i]+=temp_new[i][j];
        }
        g_ratio=g_ratio*t_new[i]/t_old[i];
    }
    g_ratio=g_ratio
           *exp(-lambda*fabs(beta_new))
           /exp(-lambda*fabs(beta[which_beta]))
           *exp(den_old-den_new);

    for(i=0;i<N;++i)
    {
        free(temp_old[i]);
        free(temp_new[i]);
    }
    free(temp_old);
    free(temp_new);

    SD = sqrt(fabs(beta_new));
    f_old = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
    SD = sqrt(fabs(beta[which_beta]));
    f_new = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);

    accept_prob = g_ratio*f_old/f_new;

    if (accept_prob > 1)
        update = 1;
    else
    {
        GetRNGstate();
        update = rbinom(1,accept_prob);
        PutRNGstate();
    }
    if (update == 1)
    {
        beta[which_beta] = beta_new;
        if (count == 1) beta_accept_count[which_beta]++;
    }
}

double update_lambda_dep(double *beta, double a, double b, int x_length)
{
    int i;
    double lambda, beta_abs[x_length+1];

    for (i=0; i<x_length+1; ++i)
        beta_abs[i] = fabs(beta[i]);
    GetRNGstate();
    lambda=rgamma((double) a+x_length+1, 1/(sum_dep(beta_abs, x_length+1)+b));
    PutRNGstate();

    return lambda;
}

void update_b_par_dep(double **b, int which_b_first, int which_b_second, double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int **b_par_accept_count, int count)
{
    int i, j, k, update, first, second, third;
    double a_old[h_length][h_length][len_E], a_new[h_length][h_length][len_E], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double b_new, b_new_vec[h_length-1][len_dummy_dep+1], min_f_old_each[len_dep], min_f_new_each[len_dep], min_f_old, min_f_new, den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;
    double **freq_new, *den_freq_new;

    GetRNGstate();

    b_new = gen_double_exp_dep(b[which_b_first][which_b_second], sqrt(fabs(b[which_b_first][which_b_second])));
    for (i=0; i<h_length-1; ++i)
    {
        for (j=0; j<len_dummy_dep+1; ++j)
        {b_new_vec[i][j]=b[i][j];}
    }
    b_new_vec[which_b_first][which_b_second]=b_new;

    freq_new=calloc(len_dep, sizeof(double*));
    for(i=0;i<len_dep;++i)
    {
        freq_new[i]=calloc(h_length, sizeof(double));
        for(j=0;j<h_length-1;++j)
        {
            freq_new[i][j]=exp(b_new_vec[j][0]);
            for(k=1;k<len_dummy_dep+1;++k)
                freq_new[i][j]*=exp(b_new_vec[j][k]*value_dep[i][k-1]);
        }
        freq_new[i][h_length-1]=1;
    }

    den_freq_new=calloc(len_dep, sizeof(double));
    for(i=0;i<len_dep;++i)
    {
        den_freq_new[i]=sum_dep(freq_new[i],h_length);
        for(j=0;j<h_length;++j)
            freq_new[i][j]/=den_freq_new[i];
        min_f_old_each[i]=find_min_dep(freq[i],h_length);
        min_f_new_each[i]=find_min_dep(freq_new[i],h_length);
    }
    min_f_old=find_min_dep(min_f_old_each, len_dep);
    min_f_new=find_min_dep(min_f_new_each, len_dep);

    assert(-min_f_old/(1-min_f_old)  < D);
    if (-min_f_new/(1-min_f_new) < D)
    {
        temp_old=calloc(N, sizeof(double*));
        temp_new=calloc(N, sizeof(double*));
        for(i=0;i<N;++i)
        {
            temp_old[i]=calloc(num_haplo_id[i], sizeof(double));
            temp_new[i]=calloc(num_haplo_id[i], sizeof(double));
        }

        g_ratio=1;
        den_old=calc_den_post_dep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_old, theta);
        den_new=calc_den_post_dep(beta, freq_new, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_new, theta);
        for(i=0;i<N;++i)
        {
            t_old[i]=0;
            t_new[i]=0;
            for(j=0;j<num_haplo_id[i];++j)
            {
                temp_old[i][j]=1;
                first=haplo_map[i][j][0]-1;
                second=haplo_map[i][j][1]-1;
                third=index_E[i]-1;
                if (y[i]==1)
                {
                    temp_old[i][j]*=theta[first][second][third];
                }
                temp_new[i][j]=temp_old[i][j]*a_new[first][second][third];
                temp_old[i][j]=temp_old[i][j]*a_old[first][second][third];
                t_old[i]+=temp_old[i][j];
                t_new[i]+=temp_new[i][j];
            }
            g_ratio=g_ratio*t_new[i]/t_old[i];

        }
        g_ratio=g_ratio*(1-min_f_new)/(1-min_f_old)*exp(den_old-den_new);
        g_ratio=g_ratio
	        *exp(-0.5*fabs(b_new))
                /exp(-0.5*fabs(b[which_b_first][which_b_second]));

        for(i=0;i<N;++i)
        {
            free(temp_old[i]);
            free(temp_new[i]);
        }
        free(temp_old);
        free(temp_new);

        SD = sqrt(fabs(b_new));
        f_old = exp(-sqrt(2)*fabs(b[which_b_first][which_b_second]-b_new)/SD)/(sqrt(2)*SD);
        SD = sqrt(fabs(b[which_b_first][which_b_second]));
        f_new = exp(-sqrt(2)*fabs(b[which_b_first][which_b_second]-b_new)/SD)/(sqrt(2)*SD);

        accept_prob = g_ratio*f_old/f_new;

        if (accept_prob > 1) update = 1;
        else update = rbinom(1, accept_prob);
        if (update == 1)
        {
            b[which_b_first][which_b_second]=b_new;
            for (i=0; i<len_dep; ++i)
            {
	        for (j=0; j<h_length; ++j)
                {
                    freq[i][j]=freq_new[i][j];
                }
            }
            if (count == 1) b_par_accept_count[which_b_first][which_b_second]++;
        }
    }

    PutRNGstate();
}

double update_D_dep(double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int len_dep, int count)
{
    int i, j, update, first, second, third;
    double a_old[h_length][h_length][len_E], a_new[h_length][h_length][len_E], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_D, min_f_each[len_dep], min_f, delta=0.05, lower, upper, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    for(i=0;i<len_dep;++i)
        min_f_each[i]=find_min_dep(freq[i],h_length);
    min_f=find_min_dep(min_f_each, len_dep);

    lower = D-delta;
    upper = D+delta;
    if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
    if (upper > 1) upper = 1;
    prop_D = runif(lower, upper);

    temp_old=calloc(N, sizeof(double*));
    temp_new=calloc(N, sizeof(double*));
    for(i=0;i<N;++i)
    {
        temp_old[i]=calloc(num_haplo_id[i], sizeof(double));
        temp_new[i]=calloc(num_haplo_id[i], sizeof(double));
    }

    g_ratio=1;
    den_old=calc_den_post_dep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_old, theta);
    den_new=calc_den_post_dep(beta, freq, prop_D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_new, theta);
    for(i=0;i<N;++i)
    {
        t_old[i]=0;
        t_new[i]=0;
        for(j=0;j<num_haplo_id[i];++j)
        {
            temp_old[i][j]=1;
            first=haplo_map[i][j][0]-1;
            second=haplo_map[i][j][1]-1;
            third=index_E[i]-1;
            if (y[i]==1)
            {
                temp_old[i][j]*=theta[first][second][third];
            }
            temp_new[i][j]=temp_old[i][j]*a_new[first][second][third];
            temp_old[i][j]=temp_old[i][j]*a_old[first][second][third];
            t_old[i]+=temp_old[i][j];
            t_new[i]+=temp_new[i][j];
        }
        g_ratio=g_ratio*t_new[i]/t_old[i];
    }
    g_ratio=g_ratio*exp(den_old-den_new);

    for(i=0;i<N;++i)
    {
        free(temp_old[i]);
        free(temp_new[i]);
    }
    free(temp_old);
    free(temp_new);

    f_new = 1/(upper-lower);
    lower = prop_D-delta;
    upper = prop_D+delta;
    if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
    if (upper > 1) upper = 1;
    f_old = 1/(upper-lower);

    accept_prob = g_ratio*f_old/f_new;
    assert(-min_f/(1-min_f)  < D);
    assert(-min_f/(1-min_f)  < prop_D);
    if (accept_prob > 1) update = 1;
    else update = rbinom(1, accept_prob);
    if (update == 1)
    {
        return prop_D;
    }
    else
        return D;

    PutRNGstate();
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, b_par, and D */
double calc_den_post_dep(double *beta, double **freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, double a[h_length][h_length][len_E], double theta[h_length][h_length][len_E])
{
    int i, j, k, l, first, second, third;
    double den=0, term[tot_uniq_mat], temp_a, temp_theta;
    for(j=0;j<len_E;++j)
    {
    for (i=0; i<tot_uniq_mat; ++i)
    {
        first=uniq_map[i][0]-1;
        second=uniq_map[i][1]-1;
        third=j;
        temp_a=calc_a_dep(freq[index_Etodep[j]-1], uniq_map[i], D);
        a[first][second][third]=temp_a;
        temp_theta=exp(beta[0]);
        for (l=0; l<len_dummy; ++l)
            temp_theta=temp_theta*exp(beta[(h_length-1)*(len_dummy+1)+l+1]*value_E[j][l]);
        if(first<h_length-1 && second<h_length-1){
            temp_theta=temp_theta*exp(beta[first+1]+beta[second+1]);
            for (k=0; k<len_dummy; ++k)
                temp_theta=temp_theta*exp((beta[first+1+(h_length-1)*(k+1)]+beta[second+1+(h_length-1)*(k+1)])*value_E[j][k]);
        }else if(first==h_length-1 && second<h_length-1){
            temp_theta=temp_theta*exp(beta[second+1]);
            for (k=0; k<len_dummy; ++k)
                temp_theta=temp_theta*exp((beta[second+1+(h_length-1)*(k+1)])*value_E[j][k]);
        }else if(first<h_length-1 && second==h_length-1){
            temp_theta=temp_theta*exp(beta[first+1]);
            for (k=0; k<len_dummy; ++k)
                temp_theta=temp_theta*exp((beta[first+1+(h_length-1)*(k+1)])*value_E[j][k]);
        }
        theta[first][second][third]=temp_theta;
        term[i]=temp_a*temp_theta;
    }
    den = den+num_E[j]*log(1+sum_dep(term, tot_uniq_mat));
    }
    return den;
}

/* Calculate a(F) that is in the denominator of the likelihood*/
double calc_a_dep(double *freq, int *per_freq, double D)
{
    int i, j, k;
    double a;
    i=per_freq[0];
    j=per_freq[1];
    if (i==j)
        {k=1;}
    else
        {k=0;}
    a = k*D*freq[i-1]+(2-k)*(1-D)*freq[i-1]*freq[j-1];
    return a;
}

/* function to find sum_dep of real numbers */
double sum_dep(double *x, int n)
{
    double sum_dep=0.0;
    int i;
    for (i=0; i<n ; ++i)
        sum_dep = sum_dep + x[i];
    return sum_dep;
}

/* function to calculate min. of an array of numbers of length n */
double find_min_dep(double *arr, int n)
{
    int i;
    double min=arr[0];
    for(i=1;i<n; ++i)
    {
        if (min > arr[i])
            min = arr[i];
    }
    return min;
}

/* function to generate from double exponential distribution */
double gen_double_exp_dep(double mean, double SD)
{
    double x, gen_exp;
    GetRNGstate();
    x = runif(0,1);
    gen_exp = rexp(SD/sqrt(2));
    PutRNGstate();
    if (x > 0.5)
        return gen_exp+mean;
    else
        return -gen_exp+mean;
}

/* function to generate from Dirichet distribution */
void dirichlet_dep(double *param, int dim, double *gen_sample)
{
    int i;
    double gen_gamma[dim], sum_dep_gamma;
    GetRNGstate();
    for (i=0; i<dim; ++i)
    {
        assert(param[i]>0);
        gen_gamma[i] = rgamma(param[i], 1);
        if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
    }
    sum_dep_gamma = sum_dep(gen_gamma, dim);
    for (i=0; i<dim; ++i)
    {
        gen_sample[i] = gen_gamma[i]/sum_dep_gamma;
    }
    PutRNGstate();
}
