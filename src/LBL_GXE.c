#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>

double calc_a(double *freq, int *per_freq, double D);

double calc_den_post1(double *beta, double *freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, double a[h_length][h_length], double theta[h_length][h_length][len_E]);

double calc_den_post2(double *beta, double **freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, double a[h_length][h_length][len_E], double theta[h_length][h_length][len_E]);

double update_lambda(double *beta, double a, double b, int x_length);

/*****************For Model 1(The main paper method)*******************/
void update_beta1(double *beta, int which_beta, double lambda, double *freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

void update_freq(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

double update_D1(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy);

/*****************For Model 2(The new method)**************************/
void update_beta2(double *beta, int which_beta, double lambda, double **freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep);

void update_b_par(double **b_par, int which_b_first, int which_b_second, double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep);

double update_D2(double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int len_dep);

/*****************From Model 1 to Model 2******************************/
void RJ12(double *beta, double *freq_indep, double **b_par, double **freq_dep, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int *model, int *skip);

/*****************From Model 2 to Model 1******************************/
void RJ21(double *beta, double *freq_indep, double **b_par, double **freq_dep, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int *model, int *skip);

double gen_double_exp(double mean, double SD);

void dirichlet(double *param, int dim, double *gen_sample);

double sum(double *x, int n);

double find_min(double *arr, int n);

void mcmc(int *y, int *N, int *n_tot, int *num_haplo_id, int *x_length, int *h_length, int *haplotype_map, int *tot_uniq_mat, int *unique_map, int *index_E, int *value_Environment, int *num_E, int *len_E, int *len_dummy, int *index_Etodep, int *value_dependent, int *len_dep, int *len_dummy_dep, double *beta, double *lambda, double *freq_indep, double *D, double *a, double *b, double *beta_out, double *lambda_out, double *freq_out, double *freq_out_indep, double *freq_out_dep, double *b_out_dep, double *D_out, int *NUM_IT, int *BURN_IN, int *count_model_1, int *count_model_2)
{

/*  Summary of notations
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
    value_dep[i][j]: all unique rows of "dep_mat" (i<*len_dep, j<*len_dummy_dep)
    *len_dep: num of unique rows of dep_mat
*/

    int i, j, k, l, m, n, ***haplo_map, h_mat[*n_tot][2], **uniq_map, **value_E, **value_dep, n_cases=0, which_beta, which_b_first, which_b_second, it1=0, it2=0, it3=0, it3_indep=0, it3_dep=0, it4_dep=0;

    int former_model;

    int *model=calloc(1, sizeof(int)), *skip=calloc(1, sizeof(int));

    double **b_par=calloc(*h_length-1, sizeof(double*));
    for(i=0;i<*h_length-1;++i)
        b_par[i]=calloc(*len_dummy_dep+1, sizeof(double));

    double **freq_dep=calloc(*len_dep, sizeof(double*));
    for(i=0;i<*len_dep;++i)
        freq_dep[i]=calloc(*h_length, sizeof(double));

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

    /*************************start MCMC here *****************************/

    /**************Do RJMCMC for 1st iteration*****************************/
    RJ12(beta, freq_indep, b_par, freq_dep, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, value_dep, *len_dep, *len_dummy_dep, model, skip);

    for (n=0; n<*NUM_IT; ++n)
    {
        if (*model == 0)
        {
            if (n >= *BURN_IN) (*count_model_1)++;

            /*****************For Model 1(The main paper method)*******************/
            if (*skip == 0){
            /* update beta parameters */
            for (i=0; i<*x_length+1; ++i)
            {
                which_beta=i;
                update_beta1(beta, which_beta, *lambda, freq_indep, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);
            }
            /* update lambda */
            *lambda = update_lambda(beta, *a, *b, *x_length);
            /* update frequencies*/
            update_freq(freq_indep, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);
            /* update D parameter */
            *D = update_D1(freq_indep, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy);}

            /*****************From Model 1 to Model 2******************************/
            former_model=0;
            RJ12(beta, freq_indep, b_par, freq_dep, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, value_dep, *len_dep, *len_dummy_dep, model, skip);
        }else{
            if (n >= *BURN_IN) (*count_model_2)++;

            /*****************For Model 2(The new method)**************************/
            if (*skip == 0){
            /* update beta parameters */
            for (i=0; i<*x_length+1; ++i)
            {
                which_beta=i;
                update_beta2(beta, which_beta, *lambda, freq_dep, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep);
            }
            /* update lambda */
            *lambda = update_lambda(beta, *a, *b, *x_length);
            /* update b and frequencies */
            for (i=0; i<*h_length-1; ++i)
            {
                for (j=0; j<*len_dummy_dep+1; ++j)
                {
                    which_b_first=i;
                    which_b_second=j;
                    update_b_par(b_par, which_b_first, which_b_second, freq_dep, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, value_dep, *len_dep, *len_dummy_dep);
                }
            }
            /* update D parameter */
            *D = update_D2(freq_dep, beta, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, *len_dep);}

            /*****************From Model 2 to Model 1******************************/
            former_model=1;
            RJ21(beta, freq_indep, b_par, freq_dep, *D, y, *N, num_haplo_id, *x_length, *h_length, haplo_map, *tot_uniq_mat, uniq_map, index_E, value_E, num_E, *len_E, *len_dummy, index_Etodep, value_dep, *len_dep, *len_dummy_dep, model, skip);
        }

        /***********output**************/
        if (n >= *BURN_IN)
        {
            for (i=0; i<*x_length+1; ++i)
            {
                beta_out[it1] = beta[i];
                ++it1;
            }

            lambda_out[it2] = *lambda;

            if (former_model == 0)
            {
                for (i=0; i<*h_length; ++i)
                {
                    freq_out_indep[it3_indep]=freq_indep[i];
                    ++it3_indep;
                    freq_out[it3]=freq_indep[i];
                    ++it3;
                }
            }
            else
            {
                for (i=0; i<*len_dep; ++i)
                {
                    for (j=0; j<*h_length; ++j)
                    {
                        freq_out_dep[it3_dep]=freq_dep[i][j];
                        ++it3_dep;
                    }
                }
                for (i=0; i<*h_length; ++i)
                {
                    freq_out[it3]=0;
                    for (j=0; j<*len_E; ++j)
                    {
                        freq_out[it3]+=freq_dep[index_Etodep[j]-1][i]*num_E[j]/(*N);
                    }
                    ++it3;
                }
                for (i=0; i<*h_length-1; ++i)
                {
                    for (j=0; j<*len_dummy_dep+1; ++j)
                    {
                        b_out_dep[it4_dep]=b_par[i][j];
                        ++it4_dep;
                    }
                }
            }

            D_out[it2] = *D;
            ++it2;
        }
    }
}

/*****************For Model 1(The main paper method)*******************/
void update_beta1(double *beta, int which_beta, double lambda, double *freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, update, first, second, third;
    double a[h_length][h_length], theta_old[h_length][h_length][len_E], theta_new[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double beta_new, beta_new_vec[x_length+1], den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;

    beta_new = gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
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
    den_old=calc_den_post1(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a, theta_old);
    den_new=calc_den_post1(beta_new_vec, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a, theta_new);
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

double update_lambda(double *beta, double a, double b, int x_length)
{
    double lambda, beta_abs[x_length+1];
    int i;

    for (i=0; i<x_length+1; ++i)
        beta_abs[i] = fabs(beta[i]);
    GetRNGstate();
    lambda=rgamma((double) a+x_length+1, 1/(sum(beta_abs, x_length+1)+b));
    PutRNGstate();

    return lambda;
}

void update_freq(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, C=10000, update, first, second, third;
    double a_old[h_length][h_length], a_new[h_length][h_length], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_freq_with_last[h_length], b_old[h_length], b_new[h_length], min_f_old, min_f_new, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    for (i=0; i<h_length; ++i)
    {
        b_old[i] = freq[i]*C;
    }
    dirichlet(b_old, h_length, prop_freq_with_last);
    min_f_old = find_min(freq, h_length);
    min_f_new = find_min(prop_freq_with_last, h_length);

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
        den_old=calc_den_post1(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_old, theta);
        den_new=calc_den_post1(beta, prop_freq_with_last, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_new, theta);
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

double update_D1(double *freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy)
{
    int i, j, update, first, second, third;
    double a_old[h_length][h_length], a_new[h_length][h_length], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_D, min_f, delta=0.05, lower, upper, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    min_f = find_min(freq, h_length);
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
    den_old=calc_den_post1(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_old, theta);
    den_new=calc_den_post1(beta, freq, prop_D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_new, theta);
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

/*****************From Model 1 to Model 2******************************/
void RJ12(double *beta, double *freq_indep, double **b, double **freq_dep, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int *model, int *skip)
{
    int i, j, k, first, second, third;
    double a_m1[h_length][h_length], a_m2[h_length][h_length][len_E], theta_m1[h_length][h_length][len_E], theta_m2[h_length][h_length][len_E], **temp_m1, **temp_m2, t_m1[N], t_m2[N];
    double min_f_m1, min_f_m2, min_f_m2_each[len_dep], den_m1, den_m2, g_ratio, accept_prob;
    double *den_freq_dep;
    *skip=0;

    GetRNGstate();

    *model = rbinom(1, 0.7);
    if (*model == 1)
    {
        for(i=0;i<h_length-1;++i)
        {
            b[i][0]=log(freq_indep[i]/freq_indep[h_length-1]);
            for(j=1;j<len_dummy_dep+1;++j)
            {
                b[i][j]=gen_double_exp(0, sqrt(0.0008));
            }
        }

        for(i=0;i<len_dep;++i)
        {
            for(j=0;j<h_length-1;++j)
            {
                freq_dep[i][j]=exp(b[j][0]);
                for(k=1;k<len_dummy_dep+1;++k)
                    freq_dep[i][j]*=exp(b[j][k]*value_dep[i][k-1]);
            }
            freq_dep[i][h_length-1]=1;
        }

        den_freq_dep=calloc(len_dep, sizeof(double));
        for(i=0;i<len_dep;++i)
        {
            den_freq_dep[i]=sum(freq_dep[i],h_length);
            for(j=0;j<h_length;++j)
                freq_dep[i][j]/=den_freq_dep[i];
            min_f_m2_each[i]=find_min(freq_dep[i],h_length);
        }
        min_f_m2=find_min(min_f_m2_each, len_dep);
        min_f_m1=find_min(freq_indep, h_length);

        temp_m1=calloc(N, sizeof(double*));
        temp_m2=calloc(N, sizeof(double*));
        for(i=0;i<N;++i)
        {
            temp_m1[i]=calloc(num_haplo_id[i], sizeof(double));
            temp_m2[i]=calloc(num_haplo_id[i], sizeof(double));
        }

        g_ratio=1;
        den_m1=calc_den_post1(beta, freq_indep, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_m1, theta_m1);
        den_m2=calc_den_post2(beta, freq_dep, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_m2, theta_m2);
        for(i=0;i<N;++i)
        {
            t_m1[i]=0;
            t_m2[i]=0;
            for(j=0;j<num_haplo_id[i];++j)
            {
                temp_m1[i][j]=1;
                temp_m2[i][j]=1;
                first=haplo_map[i][j][0]-1;
                second=haplo_map[i][j][1]-1;
                third=index_E[i]-1;
                if (y[i]==1)
                {
                    temp_m1[i][j]*=theta_m1[first][second][third];
                    temp_m2[i][j]*=theta_m2[first][second][third];
                }
                temp_m1[i][j]=temp_m1[i][j]*a_m1[first][second];
                temp_m2[i][j]=temp_m2[i][j]*a_m2[first][second][third];
                t_m1[i]+=temp_m1[i][j];
                t_m2[i]+=temp_m2[i][j];
            }
            g_ratio=g_ratio*t_m2[i]/t_m1[i];
        }
        g_ratio=g_ratio*exp(den_m1-den_m2);

        int factorial=1;
        for(i=0;i<h_length-1;++i)factorial*=i+1;

        accept_prob=g_ratio*(1-min_f_m2)/(1-min_f_m1);
        accept_prob=accept_prob/factorial/freq_indep[h_length-1]*0.3/0.7;
        for(i=0;i<h_length-1;++i)
            accept_prob*=0.25*exp(-0.5*fabs(b[i][0]))/freq_indep[i];

        for(i=0;i<N;++i)
        {
            free(temp_m1[i]);
            free(temp_m2[i]);
        }
        free(temp_m1);
        free(temp_m2);

        if (-min_f_m2/(1-min_f_m2) >= D) accept_prob = 0;
        if (accept_prob > 1) accept_prob = 1;

        *model = rbinom(1, accept_prob);
        //if (*model == 0) *skip = 1;
        //if (*model == 1) *skip = 1;
        *skip = 1;
    }

    PutRNGstate();
}

/*****************For Model 2(The new method)**************************/
void update_beta2(double *beta, int which_beta, double lambda, double **freq, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep)
{
    int i, j, update, first, second, third;
    double a[h_length][h_length][len_E], theta_old[h_length][h_length][len_E], theta_new[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double beta_new, beta_new_vec[x_length+1], den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;

    beta_new = gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
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
    den_old=calc_den_post2(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a, theta_old);
    den_new=calc_den_post2(beta_new_vec, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a, theta_new);
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
    if (update == 1) {beta[which_beta] = beta_new;}
}

void update_b_par(double **b, int which_b_first, int which_b_second, double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep)
{
    int i, j, k, update, first, second, third;
    double a_old[h_length][h_length][len_E], a_new[h_length][h_length][len_E], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double b_new, b_new_vec[h_length-1][len_dummy_dep+1], min_f_old_each[len_dep], min_f_new_each[len_dep], min_f_old, min_f_new, den_old, den_new, g_ratio, f_old, f_new, SD, accept_prob;
    double **freq_new, *den_freq_new;

    GetRNGstate();

    b_new = gen_double_exp(b[which_b_first][which_b_second], sqrt(fabs(b[which_b_first][which_b_second])));
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
        den_freq_new[i]=sum(freq_new[i],h_length);
        for(j=0;j<h_length;++j)
            freq_new[i][j]/=den_freq_new[i];
        min_f_old_each[i]=find_min(freq[i],h_length);
        min_f_new_each[i]=find_min(freq_new[i],h_length);
    }
    min_f_old=find_min(min_f_old_each, len_dep);
    min_f_new=find_min(min_f_new_each, len_dep);

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
        den_old=calc_den_post2(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_old, theta);
        den_new=calc_den_post2(beta, freq_new, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_new, theta);
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
        }
    }

    PutRNGstate();
}

double update_D2(double **freq, double *beta, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int len_dep)
{
    int i, j, update, first, second, third;
    double a_old[h_length][h_length][len_E], a_new[h_length][h_length][len_E], theta[h_length][h_length][len_E], **temp_old, **temp_new, t_old[N], t_new[N];
    double prop_D, min_f_each[len_dep], min_f, delta=0.05, lower, upper, den_old, den_new, g_ratio, f_old, f_new, accept_prob;

    GetRNGstate();

    for(i=0;i<len_dep;++i)
        min_f_each[i]=find_min(freq[i],h_length);
    min_f=find_min(min_f_each, len_dep);

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
    den_old=calc_den_post2(beta, freq, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_old, theta);
    den_new=calc_den_post2(beta, freq, prop_D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_new, theta);
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
        return prop_D;
    else
        return D;

    PutRNGstate();
}

/*****************From Model 2 to Model 1******************************/
void RJ21(double *beta, double *freq_indep, double **b, double **freq_dep, double D, int *y, int N, int *num_haplo_id, int x_length, int h_length, int ***haplo_map, int tot_uniq_mat, int **uniq_map, int *index_E, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, int **value_dep, int len_dep, int len_dummy_dep, int *model, int *skip)
{
    int i, j, first, second, third;
    double a_m1[h_length][h_length], a_m2[h_length][h_length][len_E], theta_m1[h_length][h_length][len_E], theta_m2[h_length][h_length][len_E], **temp_m1, **temp_m2, t_m1[N], t_m2[N];
    double min_f_m1, min_f_m2, min_f_m2_each[len_dep], den_m1, den_m2, g_ratio, accept_prob;
    *skip=0;

    GetRNGstate();

    *model = 1-rbinom(1, 0.3);
    if (*model == 0)
    {
        for(i=0;i<h_length;++i)
        {
            freq_indep[i]=freq_dep[0][i];
        }

        min_f_m1=find_min(freq_indep, h_length);
        for(i=0;i<len_dep;++i)
            min_f_m2_each[i]=find_min(freq_dep[i],h_length);
        min_f_m2=find_min(min_f_m2_each, len_dep);

        temp_m1=calloc(N, sizeof(double*));
        temp_m2=calloc(N, sizeof(double*));
        for(i=0;i<N;++i)
        {
            temp_m1[i]=calloc(num_haplo_id[i], sizeof(double));
            temp_m2[i]=calloc(num_haplo_id[i], sizeof(double));
        }

        g_ratio=1;
        den_m1=calc_den_post1(beta, freq_indep, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, a_m1, theta_m1);
        den_m2=calc_den_post2(beta, freq_dep, D, x_length, h_length, tot_uniq_mat, uniq_map, value_E, num_E, len_E, len_dummy, index_Etodep, a_m2, theta_m2);
        for(i=0;i<N;++i)
        {
            t_m1[i]=0;
            t_m2[i]=0;
            for(j=0;j<num_haplo_id[i];++j)
            {
                temp_m1[i][j]=1;
                temp_m2[i][j]=1;
                first=haplo_map[i][j][0]-1;
                second=haplo_map[i][j][1]-1;
                third=index_E[i]-1;
                if (y[i]==1)
                {
                    temp_m1[i][j]*=theta_m1[first][second][third];
                    temp_m2[i][j]*=theta_m2[first][second][third];
                }
                temp_m1[i][j]=temp_m1[i][j]*a_m1[first][second];
                temp_m2[i][j]=temp_m2[i][j]*a_m2[first][second][third];
                t_m1[i]+=temp_m1[i][j];
                t_m2[i]+=temp_m2[i][j];
            }
            g_ratio=g_ratio*t_m1[i]/t_m2[i];
        }
        g_ratio=g_ratio*exp(den_m2-den_m1);

        int factorial=1;
        for(i=0;i<h_length-1;++i)factorial*=i+1;

        accept_prob=g_ratio*(1-min_f_m1)/(1-min_f_m2);
        accept_prob=accept_prob*factorial*freq_indep[h_length-1]*0.7/0.3;
        for(i=0;i<h_length-1;++i)
            accept_prob/=0.25*exp(-0.5*fabs(b[i][0]))/freq_indep[i];

        for(i=0;i<N;++i)
        {
            free(temp_m1[i]);
            free(temp_m2[i]);
        }
        free(temp_m1);
        free(temp_m2);

        if (-min_f_m1/(1-min_f_m1) >= D) accept_prob = 0;
        if (accept_prob > 1) accept_prob = 1;

        *model = 1-rbinom(1, accept_prob);
        //if (*model == 1) *skip = 1;
        //if (*model == 0) *skip = 1;
        *skip = 1;
    }

    PutRNGstate();
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, frequencies, and D in Model 1*/
double calc_den_post1(double *beta, double *freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, double a[h_length][h_length], double theta[h_length][h_length][len_E])
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
        temp_a=calc_a(freq, uniq_map[i], D);
        a[first][second]=temp_a;
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
    den = den+num_E[j]*log(1+sum(term, tot_uniq_mat));
    }
    return den;
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, frequencies, and D in Model 2*/
double calc_den_post2(double *beta, double **freq, double D, int x_length, int h_length, int tot_uniq_mat, int **uniq_map, int **value_E, int *num_E, int len_E, int len_dummy, int *index_Etodep, double a[h_length][h_length][len_E], double theta[h_length][h_length][len_E])
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
        temp_a=calc_a(freq[index_Etodep[j]-1], uniq_map[i], D);
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
    den = den+num_E[j]*log(1+sum(term, tot_uniq_mat));
    }
    return den;
}

/* Calculate a(F) that is in the denominator of the likelihood*/
double calc_a(double *freq, int *per_freq, double D)
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

/* function to find sum of real numbers */
double sum(double *x, int n)
{
    double sum=0.0;
    int i;
    for (i=0; i<n ; ++i)
        sum = sum + x[i];
    return sum;
}

/* function to calculate min. of an array of numbers of length n */
double find_min(double *arr, int n)
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
double gen_double_exp(double mean, double SD)
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
void dirichlet(double *param, int dim, double *gen_sample)
{
    int i;
    double gen_gamma[dim], sum_gamma;
    GetRNGstate();
    for (i=0; i<dim; ++i)
    {
        assert(param[i] > 0);
        gen_gamma[i] = rgamma(param[i], 1);
        if (gen_gamma[i]<0.000001) gen_gamma[i] = 0.000001;
    }
    sum_gamma = sum(gen_gamma, dim);
    for (i=0; i<dim; ++i)
    {
        gen_sample[i] = gen_gamma[i]/sum_gamma;
    }
    PutRNGstate();
}
