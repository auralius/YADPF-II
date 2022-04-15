/*!
 * Auralius Manurung
 * ME - Universitas Pertamina
 * 2021
 *
 * Lotka - Volterra Fishery Problem
 *
 * Sundstrom, O., & Guzzella, L. (2009).A Generic Dynamic Programming
 * Matlab Function. 18th IEEE International Conference on Control
 * Applications, 7, 1625–1630.https://doi.org/10.1109/CCA.2009.5281131
 * 
 */

#include <iostream>
#include "CDPA.h"

/*
* This is a function call that discretize a state variable
*/
double state_fn(int n)
{
    return (static_cast<double>(n) * 0.1);
}

/*
* This is a function call that discretize an input variable
*/
double input_fn(int n)
{
    return (static_cast<double>(n) * 1);
}

/*
* Callback function for the discretized system dynamics
*/
colvec state_update_fn(colvec& X, colvec& U, double dt)
{
    colvec X_next(MAX_NSTATES, fill::zeros);
    X_next(0) = X(0) + dt * (0.02 * (X(0) - pow(X(0), 2) / 1000.0) - U(0));
    return X_next;
}

/*
* Callback function for the stage cost
*/
double stage_cost_fn(colvec& X, colvec& U, double dt)
{    
    double J = - dt * U(0);
    return J;
}

/*
* Callback function for the terminal cost
*/
double terminal_cost_fn(colvec& X, double dt)
{
    double xf = 750.0;
    double r = 1000.0;
    double J = r * pow((X(0) - xf), 2);
    return J;
}

/*
* This is the main function
*/
int main()
{
    COptVar x, u;
    CDPA dpa(1001, 200);

    dpa.set_result_dir("./result_demo2/");

    x.define(&state_fn, 10001);
    u.define(&input_fn, 11);

    dpa.define_state_update_fn(&state_update_fn);
    dpa.define_stage_cost_fn(&stage_cost_fn);
    dpa.define_terminal_cost_fn(&terminal_cost_fn);

    dpa.register_a_state_var(&x);
    dpa.register_an_input_var(&u);

    dpa.set_dt_dyn(0.01);
    dpa.solve();
}
