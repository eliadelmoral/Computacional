#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include "gsl_rng.h"

gsl_rng *tau;

using namespace std; 

int main()
{

    extern gsl_rng *tau;
    int semilla,i;
    double x,y;

    semilla=124573;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);


    for(i=0; i<=5; i++)
    {

        x=gsl_rng_uniform(tau); //entre 0 y 1
        y=gsl_rng_uniform_int(tau,20); //aleatorio entre 0 y 19
        cout << x << endl;
        cout << y << endl;
    }

    


    return 0;
}