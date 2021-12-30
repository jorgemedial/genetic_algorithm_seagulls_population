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
    Generate_EDO_Prediction(xt, pars->x0, number_of_years, pars);
    double error, max = 0;
    for(int i=0; i<number_of_years; i++){
        error = abs(xt[i]-yt[i]);
        max = error>max? error:max;
    }
    return 1/(max*max);
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

Individual* generate_individual(void){

    Individual* ind;
    Chromosome * geno;
    ODE_Parameters * pheno;
    if ((geno = (Chromosome *) malloc(sizeof(Chromosome))) == NULL ||
        (ind = (Individual* ) malloc(sizeof(Individual))) == NULL ||
        (pheno = (ODE_Parameters *) malloc(sizeof(ODE_Parameters))) == NULL){
            printf("Cannot allocate individual \n");
            return NULL;
        }
    else{    
        geno->x0 = ULONGran(21U);
        geno->phi = ULONGran(34U);
        geno->lambda = ULONGran(25U);
        geno->mu = ULONGran(25U);
        geno->sigma = ULONGran(17U);
        geno->delta = ULONGran(15U);

        ind->geno = geno;
        ind->pheno = pheno;
        genotype_to_phenotype(ind->geno, ind->pheno);
        return ind;
    }
}


