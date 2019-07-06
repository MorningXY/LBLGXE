#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>

double calc_a_indep(double *freq, int *per_freq, double D);

double calc_den_post_indep(double *beta, double *freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, double a[h_length][h_length], double theta[h_length][h_length][len_E]);

double update_lambda_indep(double *beta, double a, double b, int x_length);

void update_beta_indep(double *beta, int which_beta, double lambda, double *freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

void update_freq_indep(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

double update_D_indep(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

double gen_double_exp_indep(double mean, double SD);

void dirichlet_indep(double *param, int dim, double *gen_sample);

double sum_indep(double *x, int n);

double find_min_indep(double *arr, int n);

int no_interaction;
void mcmc_indep(int *y, int *N, int *n_tot, int *num_haplo_id, int *x_length, int *h_length, int *haplotype_map, int *tot_uniq_mat, int *unique_map, int *index_E, int *value_Environment, int *num_E, int *len_E, int *len_dummy, double *beta, double *lambda, double *freq, double *D, double *a, double *b, double *beta_out, double *lambda_out, double *freq_out, double *D_out, int *NUM_IT, int *BURN_IN, int *nointeraction)
{

/*  sum_indepmary of notations
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
*/

    int i, j, k, l, m, n, ***haplo_map, h_mat[*n_tot][2], **uniq_map, **value_E, n_cases=0, which_beta, it1=0, it2=0, it3=0, heat_it;

    no_interaction=*nointeraction;

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

    /* separating h_mat as haplo_mat */

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

    /********************** start MCMC here ***************************/

    heat_it = 0; /* tracks # of heating iterations */

    for (n=0; n<*NUM_IT; ++n)
    {
        heat_it = heat_it+1;
        if (heat_it == 101) heat_it = 1;

        /* update beta parameters */

        for (i=0; i<*x_length+1; ++i)
        {
            which_beta=i;
            update_beta_indep(beta, which_beta, *lambda, freq, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);
        }

        /* update lambda */

        *lambda = update_lambda_indep(beta, *a, *b, *x_length);

        /* update frequencies and D */

        update_freq_indep(freq, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);

        /* update D parameter */

        *D = update_D_indep(freq, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);

        if (n >= *BURN_IN)
        {
            for (i=0; i<*x_length+1; ++i)
            {
                beta_out[it1] = beta[i];
                ++it1;
            }

            lambda_out[it2] = *lambda;

            for (i=0; i<*h_length; ++i)
            {
                freq_out[it3] = freq[i];
                ++it3;
            }

            D_out[it2] = *D;
            ++it2;
        }
    }
}

void update_beta_indep(double *beta, int which_beta, double lambda, double *freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, update, first, second, third;
    double a[h_length][h_length], theta_old[h_length][h_length][len_E], theta_new[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double beta_new, beta_new_vec[x_length+1], den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;

    beta_new = gen_double_exp_indep(beta[which_beta], sqrt(fabs(beta[which_beta])));
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
    den_old=calc_den_post_indep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a, theta_old);
    den_new=calc_den_post_indep(beta_new_vec, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a, theta_new);
    for(i=0;i<N;++i)
    {
        t_old[i]=0;
        t_new[i]=0;
        for(j=0;j<num_haplo_id[i];++j)
        {
            first=haplo_map[i][j][0]-1;
            second=haplo_map[i][j][1]-1;
            third=index_E[i]-1;
            temp_old[i][j]=a[first][second];
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
    if (update == 1) {beta[which_beta] = beta_new;}
}

double update_lambda_indep(double *beta, double a, double b, int x_length)
{
    double lambda, beta_abs[x_length+1];
    int i;

    for (i=0; i<x_length+1; ++i)
        beta_abs[i] = fabs(beta[i]);
    GetRNGstate();
    lambda=rgamma((double) a+x_length+1, 1/(sum_indep(beta_abs, x_length+1)+b));
    PutRNGstate();

    return lambda;
}

void update_freq_indep(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, C=10000, update, first, second, third;
    double a_old[h_length][h_length], a_new[h_length][h_length], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_freq_with_last[h_length], b_old[h_length], b_new[h_length], min_f_old, min_f_new, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    for (i=0; i<h_length; ++i)
    {
        b_old[i] = freq[i]*C;
    }
    dirichlet_indep(b_old, h_length, prop_freq_with_last);
    min_f_old = find_min_indep(freq, h_length);
    min_f_new = find_min_indep(prop_freq_with_last, h_length);

    assert(-min_f_old/(1-min_f_old)  < D);
    if (-min_f_new/(1-min_f_new) < D)
    {
        for (i=0; i<h_length; ++i)
        {
            b_new[i] = prop_freq_with_last[i]*C;
            assert(b_new[i] > 0);
        }

        temp_old=calloc(N, sizeof(double*));
        temp_new=calloc(N, sizeof(double*));
        for(i=0;i<N;++i)
        {
            temp_old[i]=calloc(num_haplo_id[i], sizeof(double));
            temp_new[i]=calloc(num_haplo_id[i], sizeof(double));
        }

        g_ratio=1;
        den_old=calc_den_post_indep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_old, theta);
        den_new=calc_den_post_indep(beta, prop_freq_with_last, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_new, theta);
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
                temp_new[i][j]=temp_old[i][j]*a_new[first][second];
                temp_old[i][j]=temp_old[i][j]*a_old[first][second];
                t_old[i]+=temp_old[i][j];
                t_new[i]+=temp_new[i][j];
            }
            g_ratio=g_ratio*t_new[i]/t_old[i];
        }
        g_ratio=g_ratio*(1-min_f_new)/(1-min_f_old)*exp(den_old-den_new);

        for(i=0;i<N;++i)
        {
            free(temp_old[i]);
            free(temp_new[i]);
        }
        free(temp_old);
        free(temp_new);

        f_old = lgammafn(C);
        f_new = f_old;
        for (i=0; i<h_length; ++i)
        {
            f_old += (b_new[i]-1)*log(freq[i]) - lgammafn(b_new[i]);
            f_new += (b_old[i]-1)*log(prop_freq_with_last[i]) - lgammafn(b_old[i]);
        }

        accept_prob = g_ratio*exp(f_old-f_new);

        if (accept_prob > 1) update = 1;
        else update = rbinom(1, accept_prob);
        if (update == 1)
        {
            for (i=0; i<h_length; ++i)
                freq[i] = prop_freq_with_last[i];
        }
    }

    PutRNGstate();
}

double update_D_indep(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, update, first, second, third;
    double a_old[h_length][h_length], a_new[h_length][h_length], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_D, min_f, delta=0.05, lower, upper, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    min_f = find_min_indep(freq, h_length);
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
    den_old=calc_den_post_indep(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_old, theta);
    den_new=calc_den_post_indep(beta, freq, prop_D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_new, theta);
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
            temp_new[i][j]=temp_old[i][j]*a_new[first][second];
            temp_old[i][j]=temp_old[i][j]*a_old[first][second];
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
        return prop_D;
    else
        return D;

    PutRNGstate();
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, f, and d */
double calc_den_post_indep(double *beta, double *freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, double a[h_length][h_length], double theta[h_length][h_length][len_E])
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
        temp_a=calc_a_indep(freq, uniq_map[i], D);
        a[first][second]=temp_a;
        temp_theta=exp(beta[0]);
        for (l=0; l<len_dummy; ++l)
            temp_theta=temp_theta*exp(beta[(h_length-1)*(len_dummy+1)+l+1]*value_E[j][l]);
        if(first<h_length-1 && second<h_length-1){
            temp_theta=temp_theta*exp(beta[first+1]+beta[second+1]);
            if(no_interaction==0)
            {
                for (k=0; k<len_dummy; ++k)
                    temp_theta=temp_theta*exp((beta[first+1+(h_length-1)*(k+1)]+beta[second+1+(h_length-1)*(k+1)])*value_E[j][k]);
            }
        }else if(first==h_length-1 && second<h_length-1){
            temp_theta=temp_theta*exp(beta[second+1]);
            if(no_interaction==0)
            {
                for (k=0; k<len_dummy; ++k)
                    temp_theta=temp_theta*exp((beta[second+1+(h_length-1)*(k+1)])*value_E[j][k]);
            }
        }else if(first<h_length-1 && second==h_length-1){
            temp_theta=temp_theta*exp(beta[first+1]);
            if(no_interaction==0)
            {
                for (k=0; k<len_dummy; ++k)
                    temp_theta=temp_theta*exp((beta[first+1+(h_length-1)*(k+1)])*value_E[j][k]);
            }
        }
        theta[first][second][third]=temp_theta;
        term[i]=temp_a*temp_theta;
    }
    den = den+num_E[j]*log(1+sum_indep(term, tot_uniq_mat));
    }
    return den;
}

/* Calculate a(F) that is in the denominator of the likelihood*/
double calc_a_indep(double *freq, int *per_freq, double D)
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

/* function to find sum_indep of real numbers */
double sum_indep(double *x, int n)
{
    double sum_indep=0.0;
    int i;
    for (i=0; i<n ; ++i)
        sum_indep = sum_indep + x[i];
    return sum_indep;
}

/* function to calculate min. of an array of numbers of length n */
double find_min_indep(double *arr, int n)
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
double gen_double_exp_indep(double mean, double SD)
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
void dirichlet_indep(double *param, int dim, double *gen_sample)
{
    int i;
    double gen_gamma[dim], sum_indep_gamma;
    GetRNGstate();
    for (i=0; i<dim; ++i)
    {
        assert(param[i]>0);
        gen_gamma[i] = rgamma(param[i], 1);
        if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
    }
    sum_indep_gamma = sum_indep(gen_gamma, dim);
    for (i=0; i<dim; ++i)
    {
        gen_sample[i] = gen_gamma[i]/sum_indep_gamma;
    }
    PutRNGstate();
}
