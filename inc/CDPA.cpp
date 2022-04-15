#include <omp.h>
#include "CDPA.h"

CDPA::CDPA(int nhorizons, double tend)
{   
    string log_fn = dll::program_location().filename().string() + ".log";
    log::add_file_log
    (    
        log::keywords::file_name = log_fn.c_str(),
        log::keywords::rotation_size = 10 * 1024 * 1024,
        log::keywords::time_based_rotation = log::sinks::file::rotation_at_time_point(0, 0, 0),
        log::keywords::format = "[%TimeStamp%]: %Message%"
    );

    log::core::get()->set_filter(log::trivial::severity >= log::trivial::trace);
    log::add_console_log(std::cout, log::keywords::format = "[%TimeStamp%]: %Message%");
    log::add_common_attributes();

    BOOST_LOG_TRIVIAL(info) << "HELLO!";

    _state_update_fn = NULL;
    _terminal_cost_fn = NULL;
    _stage_cost_fn = NULL;

    _nhorizons = nhorizons;
    _tend = tend;
    _dt = tend / (static_cast<double>(nhorizons) - 1.0);
    _tspan = linspace(0, tend, nhorizons);

    _dt_dyn = 0.1;

    _state_vars = new COptVar[MAX_NSTATES];
    _input_vars = new COptVar[MAX_NINPUTS];

    // Fill up with 1, assuming by default there is one node for each state and input
    _NX.ones(MAX_NSTATES);
    _NU.ones(MAX_NINPUTS);

    _nstates = 0;
    _ninputs = 0;

    _state_lb.zeros(MAX_NSTATES);
    _state_ub.zeros(MAX_NSTATES);
    _input_lb.zeros(MAX_NINPUTS);
    _input_ub.zeros(MAX_NINPUTS);

    _result_dir = "./result/";
    filesystem::create_directory(_result_dir.c_str());
}

CDPA::~CDPA()
{    
    clear();

    BOOST_LOG_TRIVIAL(info) << "BYE!";
}

void CDPA::define_state_update_fn(fntype1 f)
{
    _state_update_fn = f;
}

void CDPA::clear()
{
    delete[] _state_vars;
    delete[] _input_vars;
    delete[] _U_star_matrix;
    delete[] _U;
    delete[] _X;
}

void CDPA::define_stage_cost_fn(fntype2 f)
{
    _stage_cost_fn = f;
}

void CDPA::define_terminal_cost_fn(fntype3 f)
{
    _terminal_cost_fn = f;
}

void CDPA::register_an_input_var(COptVar* optvar)
{
    _input_vars[_ninputs] = *optvar;
    _NU(_ninputs) = optvar->length();

    _NUU = static_cast<int>(prod(conv_to<vec>::from(_NU)));

    _ninputs = _ninputs + 1;

    BOOST_LOG_TRIVIAL(info) << "New input variable"
                            << ", N: " << optvar->length()
                            << ", min: " << optvar->min()
                            << ", max: " << optvar->max()
                            << ", step: " << optvar->eval_at(1) - optvar->eval_at(0);
}

void CDPA::register_a_state_var(COptVar* optvar)
{
    _state_vars[_nstates] = *optvar;
    _NX[_nstates] = optvar->length();

    _NXX = static_cast<int>(prod(conv_to<vec>::from(_NX)));    

    _nstates = _nstates + 1;

    BOOST_LOG_TRIVIAL(info) << "New state variable"
                            << ", N: " << optvar->length()
                            << ", min: " << optvar->min()
                            << ", max: " << optvar->max()
                            << ", step: " << optvar->eval_at(1) - optvar->eval_at(0);
}

void CDPA::set_result_dir(string dir)
{
    _result_dir = dir;
    if (_result_dir.back() != '/') {
        _result_dir.append("/");
    }
    filesystem::create_directory(_result_dir.c_str());
}

void CDPA::solve()
{
    openmp_setup();
    prepare();
    build_evolution_table();

    BOOST_LOG_TRIVIAL(info) << "solve() starts";
    BOOST_LOG_TRIVIAL(info) << "\tOCP time period: " << _dt;
        
    // The stage cost is a function of both the state variables and the input
    colvec X(MAX_NSTATES, fill::value(0.0));
    colvec X_next(MAX_NSTATES, fill::value(0.0));
    colvec U(MAX_NINPUTS, fill::value(0.0));

    mat Jold;

    // The terminal cost is only a function of the state variables
    BOOST_LOG_TRIVIAL(info) << "\tCalculating terminal cost...";

    s64 ind = 0;
    #pragma omp parallel for firstprivate(X, ind)
    for (int s0 = 0; s0 < _NX(0); s0++) {
        X(0) = _X[0](s0);

        for (int s1 = 0; s1 < _NX(1); s1++) {
            X(1) = _X[1](s1);

            for (int s2 = 0; s2 < _NX(2); s2++) {
                X(2) = _X[2](s2);

                for (int s3 = 0; s3 < _NX(3); s3++) {
                    X(3) = _X[3](s3);

                    for (int s4 = 0; s4 < _NX(4); s4++) {
                        X(4) = _X[4](s4);

                        for (int s5 = 0; s5 < _NX(5); s5++) {
                            X(5) = _X[5](s5);

                            for (int s6 = 0; s6 < _NX(6); s6++) {
                                X(6) = _X[6](s6);

                                for (int s7 = 0; s7 < _NX(6); s7++) {
                                    X(7) = _X[7](s7);

                                    for (int s8 = 0; s8 < _NX(8); s8++) {
                                        X(8) = _X[8](s8);

                                        for (int s9 = 0; s9 < _NX(9); s9++) {
                                            X(9) = _X[9](s9);

                                            ind = s0
                                                + s1 * _NX(0)
                                                + s2 * _NX(0) * _NX(1)
                                                + s3 * _NX(0) * _NX(1) * _NX(2)
                                                + s4 * _NX(0) * _NX(1) * _NX(2) * _NX(3)
                                                + s5 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4)
                                                + s6 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5)
                                                + s7 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6)
                                                + s8 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7)
                                                + s9 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7) * _NX(8);

                                            _J(ind) = _terminal_cost_fn(X, _dt);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << "\tCalculating stage cost...";
    for (int k = _nhorizons - 1; k > 0; k--) {
        BOOST_LOG_TRIVIAL(info) << "\t\tStage-" << k - 1;

        Jold = _J;

        #pragma omp parallel for firstprivate(X, U, X_next)
        for (int s0 = 0; s0 < _NX(0); s0++) { // --> 1st state
            X(0) = _X[0](s0);        

            for (int s1 = 0; s1 < _NX(1); s1++) { // --> 2nd state
                X(1) = _X[1](s1);               

                for (int s2 = 0; s2 < _NX(2); s2++) { // --> 3rd state
                    X(2) = _X[2](s2);
                    
                    for (int s3 = 0; s3 < _NX(3); s3++) { // --> 4th state
                        X(3) = _X[3](s3);

                        for (int s4 = 0; s4 < _NX(4); s4++) { // --> 5th state
                            X(4) = _X[4](s4);

                            for (int s5 = 0; s5 < _NX(5); s5++) { // --> 6th state
                                X(5) = _X[5](s5);

                                for (int s6 = 0; s6 < _NX(6); s6++) { // --> 7th state                                   
                                    X(6) = _X[6](s6);

                                    for (int s7 = 0; s7 < _NX(7); s7++) { // --> 8th state
                                        X(7) = _X[7](s7);

                                       for (int s8 = 0; s8 < _NX(8); s8++) {  // --> 9th state                                            
                                            X(8) = _X[8](s8);

                                            for (int s9 = 0; s9 < _NX(9); s9++) { // --> 10th state                                               
                                                X(9) = _X[9](s9);

                                                // Linear index for the current state variable combination
                                                s64 curr_ind = s0                                                       
                                                             + s1 * _NX(0)
                                                             + s2 * _NX(0) * _NX(1)
                                                             + s3 * _NX(0) * _NX(1) * _NX(2)
                                                             + s4 * _NX(0) * _NX(1) * _NX(2) * _NX(3)
                                                             + s5 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4)
                                                             + s6 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5)
                                                             + s7 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6)
                                                             + s8 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7)
                                                             + s9 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7) * _NX(8);

                                                double J_tmp = INFINITY;

                                                for (int t0 = 0; t0 < _NU(0); t0++) { // --> 1st input                                                    
                                                    U(0) = _U[0](t0);

                                                    for (int t1 = 0; t1 < _NU(1); t1++) { // --> 1st input
                                                        U(1) = _U[1](t1);

                                                        s64 input_ind = t0 + t1 * _NU(0);
                                                        s64 next_ind = _evolution_table(input_ind, curr_ind);

                                                        //double J = _stage_cost_fn(X, U, _dt) + Jold(next_ind);
                                                        double J = _stage_cost_table(input_ind, curr_ind) + Jold(next_ind);

                                                        if (J < J_tmp) {
                                                            J_tmp = J;
                                                            _descendant_matrix(curr_ind, k - 1) = next_ind;
                                                            _J(curr_ind) = J_tmp;

                                                            for (int j = 0; j < MAX_NINPUTS; j++)
                                                                _U_star_matrix[j](curr_ind, k - 1) = U(j);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } // --> 1st input
            } // --> 2nd state
        } // --> 1st state
    }

    BOOST_LOG_TRIVIAL(info) << "\tWriting the results to disks...";
    
    save_results();

    BOOST_LOG_TRIVIAL(info) << "\t\tCheck: " << _result_dir;    
    BOOST_LOG_TRIVIAL(info) << "solve() completes";
}

/*
* -------------------------------------
* Private sections
* -------------------------------------
*/ 

inline s64 CDPA::snap(double val, double min_val, double max_val, s64 N)
{
    if (max_val == min_val)  // causes division by zero!
        return 0;
        
    s64 r = (s64)floor(((val - min_val) / (max_val - min_val) * (double)(N - (s64)1))); // from 0 to N-1

    if (r > (N - (s64)1) )
        return (N - (s64)1);
    else if (r < (s64)0)
        return ((s64)0);

    return r;
}

void CDPA::openmp_setup()
{
    // Configure the OpenMP, spare one processor
    int nprocs = omp_get_num_procs();
    omp_set_num_threads(nprocs);

    BOOST_LOG_TRIVIAL(info) << "OpenMP with " << nprocs << " processors";    
}

void CDPA::prepare()
{
    // Pre-calculate the max and min values, store the results
    for (int k = 0; k < MAX_NSTATES; k++) {
        _state_ub(k) = _state_vars[k].max();
        _state_lb(k) = _state_vars[k].min();
    }

    for (int k = 0; k < MAX_NINPUTS; k++) {
        _input_ub(k) = _input_vars[k].max();
        _input_lb(k) = _input_vars[k].min();
    }

    // Pre-call the discretization function for the state variables    
    _X = new colvec[MAX_NSTATES];
    #pragma omp for 
    for (int j = 0; j < MAX_NSTATES; j++) {
        _X[j].zeros(_NX(j));
        for (int k = 0; k < _NX(j); k++) {
            _X[j](k) = _state_vars[j].eval_at(k);
        }
    }

    // Pre-call the discretization function for the input variables
    _U = new colvec[MAX_NINPUTS];
    #pragma omp for  
    for (int j = 0; j < MAX_NINPUTS; j++) {
        _U[j].zeros(_NU(j));
        for (int k = 0; k < _NU(j); k++) {
            _U[j](k) = _input_vars[j].eval_at(k);
        }
    }    

    _U_star_matrix = new mat[MAX_NINPUTS];
    for (int k = 0; k < MAX_NINPUTS; k++)
        _U_star_matrix[k] = zeros<mat>(_NXX, _nhorizons - 1);

    _J = _J.ones(_NXX) * INFINITY; // Initialize costs with a very large number        
    _descendant_matrix.zeros(_NXX, _nhorizons);
}

void CDPA::build_evolution_table()
{
    BOOST_LOG_TRIVIAL(info) << "build_evolution_table() starts";
    BOOST_LOG_TRIVIAL(info) << "\tSystem dynamic's time period: " << _dt_dyn;

    _evolution_table.zeros(_NUU, _NXX);
    _stage_cost_table.zeros(_NUU, _NXX);

    // The stage cost is a function of both the state variables and the input
    colvec X(MAX_NSTATES, fill::value(0));
    colvec X_next(MAX_NSTATES, fill::value(0));
    colvec U(MAX_NINPUTS, fill::value(0));

    int r = static_cast<int>(_dt / _dt_dyn);
    
    #pragma omp parallel for firstprivate(X, U, X_next)
    // --> 1st state
    for (int s0 = 0; s0 < _NX(0); s0++) {
        X(0) = _X[0](s0);

        for (int s1 = 0; s1 < _NX(1); s1++) { // --> 2nd state
            X(1) = _X[1](s1);

            for (int s2 = 0; s2 < _NX(2); s2++) { // --> 3rd state
                X(2) = _X[2](s2);

                for (int s3 = 0; s3 < _NX(3); s3++) { // --> 4th state
                    X(3) = _X[3](s3);

                    for (int s4 = 0; s4 < _NX(4); s4++) { // --> 5th state
                        X(4) = _X[4](s4);

                        for (int s5 = 0; s5 < _NX(5); s5++) { // --> 6th state
                            X(5) = _X[5](s5);

                            for (int s6 = 0; s6 < _NX(6); s6++) { // --> 7th state
                                X(6) = _X[6](s6);

                                for (int s7 = 0; s7 < _NX(7); s7++) { // --> 8th state
                                    X(7) = _X[7](s7);

                                    for (int s8 = 0; s8 < _NX(8); s8++) { // --> 9th state
                                        X(8) = _X[8](s8);

                                        for (int s9 = 0; s9 < _NX(9); s9++) { // --> 10th state
                                            X(9) = _X[9](s9);

                                            // Linear index for the current state variable combination
                                            s64 state_ind = s0
                                                + s1 * _NX(0)
                                                + s2 * _NX(0) * _NX(1)
                                                + s3 * _NX(0) * _NX(1) * _NX(2)
                                                + s4 * _NX(0) * _NX(1) * _NX(2) * _NX(3)
                                                + s5 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4)
                                                + s6 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5)
                                                + s7 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6)
                                                + s8 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7)
                                                + s9 * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7) * _NX(8);

                                            for (int t0 = 0; t0 < _NU(0); t0++) { // --> 1st input
                                                U(0) = _U[0](t0);

                                                for (int t1 = 0; t1 < _NU(1); t1++) { // --> 1st input
                                                    U(1) = _U[1](t1);

                                                    s64 input_ind = t0 + t1 * _NU(0);                                                      
                                                    
                                                    X_next = X;
                                                    for (int tn = 0; tn < r; tn++)
                                                        X_next = _state_update_fn(X_next, U, _dt_dyn);

                                                    // Make sure they are snapped to the defined grids and are inside the boundaries        
                                                    s64 n[MAX_NSTATES];
                                                    for (int j = 0; j < MAX_NSTATES; j++)
                                                        n[j] = snap(X_next(j), _state_lb(j), _state_ub(j), _NX(j));

                                                    // Find the corresponding linear index
                                                    s64 nextstate_ind = n[0] +
                                                                        n[1] * _NX(0) +
                                                                        n[2] * _NX(0) * _NX(1) +
                                                                        n[3] * _NX(0) * _NX(1) * _NX(2) +
                                                                        n[4] * _NX(0) * _NX(1) * _NX(2) * _NX(3) +
                                                                        n[5] * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) +
                                                                        n[6] * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) +
                                                                        n[7] * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) +
                                                                        n[8] * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7) +
                                                                        n[9] * _NX(0) * _NX(1) * _NX(2) * _NX(3) * _NX(4) * _NX(5) * _NX(6) * _NX(7) * _NX(8);

                                                    _evolution_table(input_ind, state_ind) = nextstate_ind;
                                                    _stage_cost_table(input_ind, state_ind) = _stage_cost_fn(X, U, _dt);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } // --> 1st input
        } // --> 2nd state
    } // --> 1st state
    
    BOOST_LOG_TRIVIAL(info) << "build_evolution_table() completes";
}

void CDPA::save_results()
{
    string target;
    // The descendet matrix, we will trace this matrix with MATLAB / OCTAVE
    target = _result_dir + "descendent_matrix.bin";
    _descendant_matrix.save(target.c_str(), raw_binary);

    // The optimal policies, we will trace this matrix with MATLAB / OCTAVE
    target = _result_dir + "u0_star_matrix.bin";
    _U_star_matrix[0].save(target.c_str(), raw_binary);

    target = _result_dir + "u1_star_matrix.bin";
    _U_star_matrix[1].save(target.c_str(), raw_binary);

    // Time
    target = _result_dir + "tspan.bin";
    _tspan.save(target.c_str(), raw_binary);

    // Number of nodes for each state
    target = _result_dir + "NX.bin";
    _NX.save(target.c_str(), raw_binary);

    // Number of nodes for each input
    target = _result_dir + "NU.bin";
    _NU.save(target.c_str(), raw_binary);

    // Upper boundary of each state
    target = _result_dir + "state_ub.bin";
    _state_ub.save(target.c_str(), raw_binary);

    // Lower boundary of each input
    target = _result_dir + "input_lb.bin";
    _input_lb.save(target.c_str(), raw_binary);

    // Upper boundary of each input
    target = _result_dir + "input_ub.bin";
    _input_ub.save(target.c_str(), raw_binary);

    // Lower boundary of each state
    target = _result_dir + "state_lb.bin";
    _state_lb.save(target.c_str(), raw_binary);

    // The discretized state variables 
    target = _result_dir + "X0.bin";
    _X[0].save(target.c_str(), raw_binary);

    target = _result_dir + "X1.bin";
    _X[1].save(target.c_str(), raw_binary);

    target = _result_dir + "X2.bin";
    _X[2].save(target.c_str(), raw_binary);

    target = _result_dir + "X3.bin";
    _X[3].save(target.c_str(), raw_binary);

    target = _result_dir + "X4.bin";
    _X[4].save(target.c_str(), raw_binary);

    target = _result_dir + "X5.bin";
    _X[5].save(target.c_str(), raw_binary);

    target = _result_dir + "X6.bin";
    _X[6].save(target.c_str(), raw_binary);

    target = _result_dir + "X7.bin";
    _X[7].save(target.c_str(), raw_binary);

    target = _result_dir + "X8.bin";
    _X[8].save(target.c_str(), raw_binary);

    target = _result_dir + "X9.bin";
    _X[9].save(target.c_str(), raw_binary);

    // The discretized input variables
    target = _result_dir + "U0.bin";
    _U[0].save(target.c_str(), raw_binary);

    target = _result_dir + "U1.bin";
    _U[1].save(target.c_str(), raw_binary);

    target = _result_dir + "J.bin";
    _J.save(target.c_str(), raw_binary);
}