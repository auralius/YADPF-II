/**
 * \file CDPA.h
 * \brief The main class for the implementation of dynamic programming algorithm.
 */

#pragma once

#define ARMA_USE_LAPACK

#include <iostream>
#include <fstream>
#include <utility>
#include <string>

#include <boost/filesystem.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#include <boost/dll.hpp>

#include <armadillo>

#include "COptVar.h"

using namespace std;
using namespace arma;
using namespace boost;


typedef colvec  (*fntype1)  (colvec&, colvec&, double);  /// State update function 
typedef double  (*fntype2)  (colvec&, colvec&, double);  /// Stage cost function 
typedef double  (*fntype3)  (colvec&, double);           /// Terminal cost funtion 

/** 
 * Definitions for maximum number of stateand input variables
 * Do not change these two definitions! 
 */
constexpr int MAX_NSTATES = 10;
constexpr int MAX_NINPUTS = 2;

/**
 * A class for the  implementation of Dynamic Programming Algorithm 
 */

class CDPA
{
public:
    CDPA(int nhorizons, double tend);
    ~CDPA();
    void clear();

    /**
     * Set the callback function for the state update process. 
     * During the state update process, we simulate the system dynamic one time step ahead. 
     * 
     * Prototype: colvec *foo(colvec &X, colvec &U, double dt);
     * 
     * @param f the address of the callback function 
     */
    void define_state_update_fn(fntype1 f);

    /**
    * Set the callback function for the stage cost function. 
     * The terminal cost function is a function of both the state and input varibles. 
     * 
     * Prototype: double foo(colvec &X, colvec &U, double dt);
     *
     * @param f the address of the callback function
     */
    void define_stage_cost_fn(fntype2 f);

    /**
     * Set the callback function for the terminal cost function. 
     * The terminal cost function is only a function of the state varibles. 
     * 
     * Prototype: double foo(colvec &X, double dt);
     *
     * @param f the address of the callback function
     */
    void define_terminal_cost_fn(fntype3 f);

    /**
     * Register an input variable.
     *
     * @param optvar the pointer to the optimization variable class
     */
    void register_an_input_var(COptVar* optvar);

    /**
     * Register a state variable.
     *
     * @param optvar the pointer to the optimization variable class
     */
    void register_a_state_var(COptVar* optvar);

    /**
     * Define the folder where the results will be stored. 
     * The root directory is: "bin".
     *
     * @param dir the folder name
     */
    void set_result_dir(string dir);

    void set_dt_dyn(double t) { _dt_dyn = t; };

    /**
     * Run the solver.
     */
    void solve();
    
private:
    /**
     * Save te results.
     */
    void save_results();

    /**
     * Configure the OpenMP.
     */
    void openmp_setup();

    /**
     * Preparation before running the solve function.
     */
    void prepare();
    
    /** 
     * Convert val to its closest integer (snap to integer). 
     * After that, the integer number is scaled from 0 to N-1 where min_val is be converted to 0 and max_val is converted to N-1. 
     * @param val value to snap
     * @param min_val the minimum value
     * @param max_val the maximum value
     * @param N the largest integer number
     * 
     * @return the closest integer to val
     */
    s64 snap(double val, double min_val, double max_val, s64 N);

    /**
     * Get number of state variables.
     * 
     * @return _nstates the number of state variables.
     */
    inline int nstates() { return _nstates; };

    /**
     * Get number of input variables
     * 
     * @return _ninputs the number of input variables
     */
    inline int ninputs() { return _ninputs; };

    /**
     * Get the horizon length. 
     * Horizon length is given by the data length of the discretized time. 
     * As an example: 0<=t<=10, discretized at dt=1, gives us a horizon length of 11.
     * 
     * @return _nhorizons the horizon length.
     */
    inline int nhorizons() { return _nhorizons; };


    /*!
     * Get the end time. The start time is alwayas 0.
     * 
     * @return _ten the end time
     */
    inline double tend() const { return _tend; };    

    /*!
     * To speed up the computation, instead of computing the system dynamics repetitively in every stage,
     * we will store the results in a matrix. Thus, we only need to compute the system dynamics once.
     */
    void build_evolution_table();

    fntype1  _state_update_fn;   /// Pointer to the state-update callback function 
    fntype2  _stage_cost_fn;     /// Pointer to the stage-cost callback function 
    fntype3  _terminal_cost_fn;  /// Pointer to the terminal-cost callback function 

    int _nstates;    /// Number of state variables
    int _ninputs;    /// Number of input variables
    int _nhorizons;  /// Horizon length

    COptVar* _input_vars;  /// Array that holds all input variables
    COptVar* _state_vars;  /// Array that holds all state variables

    s64_colvec _NX;  /// Number of nodes of each state variable
    s64_colvec _NU;  /// Number of nodes of each input variable

    s64 _NXX;  /// Total node combinations of all state variables
    s64 _NUU;  /// Total node combinations of all input variables

    double _tend;   /// End time
    double _dt;     /// OCP time period
    colvec _tspan;  /// Time vector

    double _dt_dyn; /// System dynamic's time period

    colvec _J;
    mat* _U_star_matrix;        /// Contains the optimal inputs over the time horizon
    s64_mat _descendant_matrix; /// Contains indices of the optimal states over the time horizon

    colvec _state_ub; /// Upper boundaries of each state variable
    colvec _state_lb; /// Lower boundaries of each state variable
    colvec _input_ub; /// Upper boundaries of each input variable
    colvec _input_lb; /// Lower boundaries of each input variable

    /**
     * Saome variables can be pre-calculated and then stored in the memory to speed up the iteration. 
     */
    colvec *_X; /// Holds the discretized state variables. 
    colvec *_U; /// Holds the discretized input variables.
    s64_mat _evolution_table; /// Holds the next state of the system for all possible inputs
    mat _stage_cost_table;

    string _result_dir;  /// Directory where the results will be stored    
};