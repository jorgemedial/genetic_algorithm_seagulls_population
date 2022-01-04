#define ElliotSigmoidSCALE 1000
#define TwoElliotSigmoidSCALE 2000
#include <stdio.h>
#include "RKF78.h"
#include <math.h>
#include "random.c" 



void read_data(double *yt, unsigned short number_of_years){
    FILE *file;

    file = fopen("./data/andouin_population.txt", "r");   

    short int year;
    int scanned;
    for(int i = 0; i<number_of_years && !feof(file); i++){
        scanned = fscanf(file, "%hd %lf", &year, yt+i);
        if (scanned==0) break;
    }

    fclose(file);

}

void save_data(double *yt, double *xt, unsigned short number_of_years){
    FILE *file;

    file = fopen("./results/andouin_population_prediction.txt", "w");   

    int writen;
    for(int i = 0; i<number_of_years && !feof(file); i++){
        writen = fprintf(file, "%hd %lf %lf\n", 2018-number_of_years+i, yt[i], xt[i]);
        if (writen==0) break;
    }

    fclose(file);

}

double ElliotSigmoid(double x, double sigma, double delta) {
x = sigma*(x-delta);
return x/(ElliotSigmoidSCALE + fabs(x));
}

double Psi(double x, double mu, double sigma, double delta){
    double ZeRoParsThreshold = 1e-10;
    if(fabs(sigma) < ZeRoParsThreshold) return 1.0;
    double ES = ElliotSigmoid(x, sigma, delta);
    sigma *= delta; x /= delta;
    if(x < delta) {
    ES = ES * (x + (mu*(1.0-x)*(sigma + ElliotSigmoidSCALE)) / (sigma + TwoElliotSigmoidSCALE));
    }
    return ((1 - ES)*(sigma + TwoElliotSigmoidSCALE)) / (sigma*(1+mu) + TwoElliotSigmoidSCALE);
}

typedef struct {
    double x0;
    double phi;
    double beta;
    double lambda;
    double mu;
    double sigma;
    double delta;
} ODE_Parameters;


typedef struct {
    unsigned long x0;
    unsigned long phi;
    unsigned long beta;
    unsigned long lambda;
    unsigned long mu;
    unsigned long sigma;
    unsigned long delta;
} Chromosome;

typedef struct {
    Chromosome* geno;
    ODE_Parameters* pheno;
    double fitness;
    double accum_fitness;
} Individual;

void MigrationODE(double t, double x, double *der, void *Params){
ODE_Parameters *par = (ODE_Parameters *) Params; // Pointer cast to save typing and thinking
*der = par->phi * x - par->beta*x*x - par->lambda*Psi(x, par->mu, par->sigma, par->delta);
}

#define HMAX 1.0
#define HMIN 1.e-6
#define RKTOL 1.e-8
int Generate_EDO_Prediction( double *xt, double x0,
unsigned short number_of_years,
ODE_Parameters *pars ){
    register unsigned ty;
    xt[0] = x0; // Storing IC x(0)
    for(ty=1; ty < number_of_years; ty++) xt[ty] = 0.0;
    double t = 0.0, err, h = 1.e-3;
    for(ty=1; ty < number_of_years; ty++) { 
        int status;
        while(t+h < ty) {
            status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
            if(status) return status;
            } // Adaptative stepsize h. To assure stopping at t = ty
        h = ty - t;
        status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
        if(status) return status;
        xt[ty] = x0;
        }
return 0;
}


double fitness(double *yt, double* xt, 
               double x0, unsigned short number_of_years,
               ODE_Parameters *pars){
    Generate_EDO_Prediction(xt, x0, number_of_years, pars);
    double error, max = 0.0;
    for(int i=0; i<number_of_years; i++){
        error = fabs(xt[i]-yt[i]);
        //printf("yt: %e, xt: %e, error: %e\n", yt[i], xt[i], error);
        max = max>error?max:error;
    }
    return 1/(max*max*max);
}


void genotype_to_phenotype(Chromosome * geno, ODE_Parameters * pheno){
    pheno->x0 = geno->x0*0.0079155006005766;
    pheno->phi = geno->phi*5.841138773e-9 - 100;
    pheno->lambda = geno->lambda*8.940696982762e-5;
    pheno->mu = geno->mu*5.960464655174e-7;
    pheno->sigma = geno->sigma*0.007629452739355006;
    pheno->delta = geno->delta*0.7629627368999298;
    pheno->beta = 0.000024382635446;
}

Individual* generate_individual(int empty){

    Individual* ind;
    Chromosome * geno;
    ODE_Parameters * pheno;
    if ((geno = (Chromosome *) malloc(sizeof(Chromosome))) == NULL ||
        (ind = (Individual* ) malloc(sizeof(Individual))) == NULL ||
        (pheno = (ODE_Parameters *) malloc(sizeof(ODE_Parameters))) == NULL){
            printf("Cannot allocate individual \n");
            return NULL;
        }
    if(!empty){    
        geno->x0 = ULONGran(21U);
        geno->phi = ULONGran(34U);
        geno->lambda = ULONGran(25U);
        geno->mu = ULONGran(25U);
        geno->sigma = ULONGran(17U);
        geno->delta = ULONGran(15U);
    }
    ind->geno = geno;
    ind->pheno = pheno;
    genotype_to_phenotype(ind->geno, ind->pheno);
    return ind;
}

Individual** generate_population(int n_inds, int empty){
    Individual **inds = (Individual **) malloc(sizeof(Individual)*n_inds);
    for(int i = 0; i<n_inds; i++){
        inds[i] = generate_individual(empty);
    }
    return inds;
}



void UniformCrossover(unsigned long p1, unsigned long p2,
    unsigned long *f1, unsigned long *f2, unsigned char len, float prob){
    unsigned long mask = 0U;
    register unsigned char i;
    for(i=0; i < len; i++) if(uniform() < prob) mask = mask | (1U << i);
    *f1 = (p1 & mask) | (p2 & ~mask);
    *f2 = (p2 & mask) | (p1 & ~mask);
}

void BitFlipMutation(unsigned long *f, unsigned char len, float prob){
    register unsigned char i;
    for(i=0; i < len; i++) if(uniform() < prob) *f = (*f)^(1U << i);
}

void mate(Chromosome *p1, Chromosome *p2,
    Chromosome *f1, Chromosome *f2){
        UniformCrossover(p1->x0, p2->x0, &(f1->x0), &(f2->x0), 21, 0.5);
        UniformCrossover(p1->phi, p2->phi, &(f1->phi), &(f2->phi), 34, 0.5);
        UniformCrossover(p1->lambda, p2->lambda, &(f1->lambda), &(f2->lambda), 25, 0.5);
        UniformCrossover(p1->mu, p2->mu, &(f1->mu), &(f2->mu), 25, 0.5);
        UniformCrossover(p1->sigma, p2->sigma, &(f1->sigma), &(f2->sigma), 17, 0.5);
        UniformCrossover(p1->delta, p2->delta, &(f1->delta), &(f2->delta), 15, 0.5);
}

void mutate(Chromosome *f1, float prob){
        BitFlipMutation(&(f1->x0), 21, prob);
        BitFlipMutation(&(f1->phi), 34, prob);
        BitFlipMutation(&(f1->lambda), 25, prob);
        BitFlipMutation(&(f1->mu), 25, prob);
        BitFlipMutation(&(f1->sigma), 17, prob);
        BitFlipMutation(&(f1->delta), 15, prob);
}

int roulette_selection(Individual **inds, int n_ind, double total_fitness){
    double roulette_selection = uniform()*total_fitness;
    int found = 0;
    int pos = floor(n_ind/2);
    int left_margin = 0;
    int right_margin = n_ind-1;
    while(!found){
        // printf("pos: %d, left_margin: %d, right_margin: %d\n", pos, left_margin, right_margin);
        // printf("indf: %e, roulette:%e\n", inds[pos]->accum_fitness, roulette_selection);

        if(inds[pos]->accum_fitness > roulette_selection){
            right_margin = pos>0? pos-1:0; 
        } 
        else{
            if(pos==n_ind-1){
                found = 1;
                break;
            }
            if(inds[pos+1]->accum_fitness > roulette_selection ){
                found = 1;
                break;
            }
            else left_margin = pos<n_ind-1? pos+1:n_ind-1;
        } 
        pos = (right_margin + left_margin)/2;
    }
    return pos;

}

double get_mating_chances(Individual **inds, int n_ind){
    // Compute probabilities
    double total_fitness = 0;
    for(int i=0; i<n_ind; i++){
        inds[i]->accum_fitness = total_fitness;
        total_fitness += inds[i]->fitness;
    }
    printf("average fitness: %1.9e\n", total_fitness/n_ind);
    return total_fitness;
}


int *select_parents(Individual **inds, int n_ind, double total_fitness){
    int *parents = (int*) malloc(2*sizeof(int));
    if(parents==NULL) exit(66);
    parents[0] = roulette_selection(inds, n_ind, total_fitness);
    do{
        parents[1] = roulette_selection(inds, n_ind, total_fitness);
    }while(parents[0]==parents[1]);

    return parents;
}

void update_fitness(Individual** inds, int n_inds, unsigned short number_of_years,
                    double *yt, double *xt){
    for(int i = 0; i<n_inds; i++){
        genotype_to_phenotype(inds[i]->geno, inds[i]->pheno);
        inds[i]->fitness = fitness(yt, xt, inds[i]->pheno->x0, number_of_years, inds[i]->pheno);
    }
}

void create_new_gen(Individual** inds, Individual** children, int n_inds){
    double total_fitness = get_mating_chances(inds, n_inds);
    int *parents;
    for(int i=0; i<n_inds-1; i+=2){
        parents = select_parents(inds, n_inds, total_fitness);
        // printf("i: %d, p1 index: %d, p2 index: %d\n", i, parents[0], parents[1]);
        // printf("p1: %lu, p2: %lu, c1: %lu, c2: %lu\n", inds[parents[0]]->geno->delta, inds[parents[1]]->geno->delta, 
        //                                            children[i]->geno->delta, children[i+1]->geno->delta);

        mate(inds[parents[0]]->geno, inds[parents[1]]->geno, children[i]->geno, children[i+1]->geno);
        // printf("p1: %lu, p2: %lu, c1: %lu, c2: %lu\n", inds[parents[0]]->geno->delta, inds[parents[1]]->geno->delta, 
        //                                            children[i]->geno->delta, children[i+1]->geno->delta);

        mutate(children[i]->geno, 0.01);
        mutate(children[i+1]->geno, 0.01);
        free(parents);
    }

}