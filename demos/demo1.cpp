/*!
 * Auralius Manurung
 * ME - Universitas Pertamina
 * 2022
 *
 * A mass(1 kg) is moving from x = 0 to x = 0.5 in exactly 1 second because of
 * an exteral force.The damping coefficient is 0.1.At the destination, the
 * mass must stop moving.The external force is bounded(-4 N to 4 N).
 * 
 */

#include <iostream>
#include "CDPA.h"

/*
* This is a function call that discretize a state variable
*/
double state_fn(int n)
{
    return (static_cast<double>(n) * 0.001);
}

/*
* This is a function call that discretize an input variable
*/
double input_fn(int n)
{
    return (static_cast<double>(n) *  0.1 - 4.0);
}

/*
* Callback function for the discretized system dynamics
*/
colvec state_update_fn(colvec& X, colvec& U, double dt)
{
    double m = 1.0;
    double b = 0.1;

    colvec X_next(MAX_NSTATES, fill::zeros);
    X_next(0) = X(0) + dt * X(1);
    X_next(1) = X(1) - b / m * X(1) * dt + dt / m * U(0);

    return X_next;
}

/*
* Callback function for the stage cost
*/
double stage_cost_fn(colvec&X, colvec&U, double dt)
{
    double r = 1.0;
    double J = r * pow(U(0), 2) * dt;
    return J;
}

/*
* Callback function for the terminal cost
*/
double terminal_cost_fn(colvec& X, double dt)
{
    double r2 = 1000.0;
    double r3 = 1000.0;

    double x1f = 0.5;
    double x2f = 0.0;

    double J = r2 * pow((X(0) - x1f),2) + r3 * pow((X(1) - x2f),2);
    return J;
}

/*
* This is the main function
*/
int main()
{
    COptVar x1, x2, u;
    CDPA dpa(11, 1);

    dpa.set_result_dir("./result_demo1/");

    x1.define(&state_fn, 1001);
    x2.define(&state_fn, 1001);
    u.define(&input_fn, 81);

    dpa.define_state_update_fn(&state_update_fn);
    dpa.define_stage_cost_fn(&stage_cost_fn);
    dpa.define_terminal_cost_fn(&terminal_cost_fn);

    dpa.register_a_state_var(&x1);
    dpa.register_a_state_var(&x2);
    dpa.register_an_input_var(&u);

    dpa.set_dt_dyn(0.01);
    dpa.solve();
}
