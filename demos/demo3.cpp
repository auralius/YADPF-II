/**
 * Auralius Manurung(manurung.auralius@gmail.com)
 * ME - Universitas Pertamina
 * November 2021
 *
 * Stabilization of an F8 aircraft
 *
 * Banks, S.P., & Mhana, K.J. (1992).Optimal controland stabilization
 * for nonlinear systems.IMA Journal of Mathematical Controland
 * Information, 9(2), 179–196.https://doi.org/10.1093/imamci/9.2.179
 *
 * Kaya, C.Y., & Noakes, J.L. (2003).Computational Method for
 * Time - Optimal Switching Control.Journal of Optimization Theory and
 * Applications, 117(1), 69–92.https ://doi.org/10.1023/A:1023600422807
 * 
 * Garrard, W.L., & Jordan, J.M. (1977).Design of nonlinear automatic
 * flight control systems.Automatica, 13(5), 497–505.
 * https ://doi.org/10.1016/0005-1098(77)90070-X
 *
 * Kaya, C.Y., & Noakes, J.L. (1996).Computationsand time - optimal
 * controls.Optimal Control Applications and Methods, 17(3), 171–185.
 * https ://doi.org/10.1002/(SICI)1099-1514(199607/09)17:3<171::AID-OCA571>3.0.CO;2-9
 */
#define _USE_MATH_DEFINES 

#include <math.h>
#include "CDPA.h"

#include <iostream>


inline double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

// X0
double state0_fn(int n)
{
    return (static_cast<double>(n) * 0.002 - 0.1);
}

// X1
double state1_fn(int n)
{
    return (static_cast<double>(n) * 0.002 - 0.6);
}

// X2
double state2_fn(int n)
{
    return (static_cast<double>(n) * 0.002 - 0.7);
}

// U0
double input_fn(int n)
{
    return (deg2rad(static_cast<double>(n) * 0.1 - 3.0));
    //return (deg2rad(static_cast<double>(n) *  0.2 - 3.0));
}

// Callback function for the discretized system dynamics
colvec state_update_fn(colvec& X, colvec& U, double dt)
{
    colvec Xnext(MAX_NSTATES, fill::zeros);
    
    Xnext(0) = (-0.877 * X(0) + X(2) - 0.088 * X(0) * X(2) + 0.47 * pow(X(0), 2) -
        0.019 * pow(X(1), 2) - pow(X(0), 2) * X(2) + 3.846 * pow(X(0), 3) -
        0.215 * U(0) + 0.28 * pow(X(0), 2) * U(0) + 0.47 * X(0) * pow(U(0), 2) +
        0.63 * pow(U(0), 3)) * dt + X(0);
    Xnext(1) = X(2) * dt + X(1);
    Xnext(2) = (-4.208 * X(0) - 0.396 * X(2) - 0.47 * pow(X(0), 2) - 3.564 * pow(X(0), 3) -
        20.967 * U(0) + 6.265 * pow(X(0), 2) * U(0) + 46.0 * X(0) * pow(U(0), 2) +
        61.4 * pow(U(0), 3)) * dt + X(2);

    return Xnext;
}

// Callback function for the stage cost
double stage_cost_fn(colvec &X, colvec&U, double dt)
{
    double x_bar = 0.01;
    double J = (pow(X(0)/x_bar, 4) + pow(X(1)/x_bar, 4) + pow(X(2)/x_bar, 4))/21.0;

    return J;
}

// Callback function for the terminal cost
double terminal_cost_fn(colvec& X, double dt)
{
    double x_bar = 0.01;
    double J = (pow(X(0)/x_bar, 4) + pow(X(1)/x_bar, 4) + pow(X(2)/x_bar, 4))/21.0;

    return J;
}

/*
* This is the main function
*/
int main()
{
    COptVar x0, x1, x2, u;
    CDPA dpa( 21, 10);

    dpa.set_result_dir("./result_demo3/");

    x0.define(&state0_fn, 351);
    x1.define(&state1_fn, 401);
    x2.define(&state2_fn, 551);
    u.define(&input_fn, 61);

    dpa.define_state_update_fn(&state_update_fn);
    dpa.define_stage_cost_fn(&stage_cost_fn);
    dpa.define_terminal_cost_fn(&terminal_cost_fn);

    dpa.register_a_state_var(&x0);
    dpa.register_a_state_var(&x1);
    dpa.register_a_state_var(&x2);
    dpa.register_an_input_var(&u);

    dpa.set_dt_dyn(0.01);
    dpa.solve();
}
