#include "lib/genetic.c"

int main(){

    unsigned short number_of_years = 37;
    double *xt = (double* ) malloc(number_of_years*sizeof(double));
    double *yt = (double* ) malloc(number_of_years*sizeof(double));

    read_data(yt, number_of_years);
    int n_individuals = 5000;
    Individual *inds[n_individuals];
    double error;
    for(int i = 0; i<n_individuals; i++){
        inds[i] = generate_individual();
        error = fitness(yt, xt, inds[i]->pheno->x0, number_of_years, inds[i]->pheno);
        if (error>1e-7) printf("error: %e\n", error);
    }
     

    return 0;
}