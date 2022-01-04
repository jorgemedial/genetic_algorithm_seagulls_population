#include "lib/genetic.c"

int main(){

    unsigned short number_of_years = 12;
    double *xt = (double* ) malloc(number_of_years*sizeof(double));
    double *yt = (double* ) malloc(number_of_years*sizeof(double));

    read_data(yt, number_of_years);


    ODE_Parameters test, test2;
    test.x0 = 16081.26;
    test.beta = 0.000024382635446;
    test.phi = 0.263109;
    test.lambda = 2624.9999;
    test.mu = 19.375;
    test.sigma = 7.81;
    test.delta = 18767;

    test2.x0 = 16243.62;
    test2.beta = 0.000024382635446;
    test2.phi = 0.16625;
    test2.lambda = 1500;
    test2.mu = 0.494602;
    test2.sigma = 31.25;
    test2.delta = 8593;

    Individual *random = generate_individual(0);

    int n_inds = 500;
    Individual **inds = generate_population(n_inds, 0);
    Individual **children = generate_population(n_inds, 1);
 
    
    // for(int i=0; i<number_of_years; i++){
    //     printf("yt: %e\n", yt[i]);
    // }
    // return 0;
    for(int i = 0; i<100; i++){
        update_fitness(inds, n_inds, number_of_years, yt, xt);
        create_new_gen(inds, children, n_inds);
        update_fitness(children, n_inds, number_of_years, yt, xt);
        create_new_gen(children, inds, n_inds);

    }   

    save_data(yt, xt, number_of_years);
    // for(int i = 0; i<n_inds; i++){
    //     printf("elem %d fitness: %e\n", i, inds[i]->fitness);
    // }

    return 0;
}