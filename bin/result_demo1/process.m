% -------------------------------------------------------------------------
% Trace the results from the C++ program
% This is mainly written for OCTAVE and might need a little tuning for 
% MATLAB
%
% Auralius Manurung
% ME, Universitas Pertamina, 2022
% auralius.manurung@gmail.com
% -------------------------------------------------------------------------

clc
close all
clear

% -------------------------------------------------------------------------
% Do not modify these two lines

MAX_NSTATES = 10;
MAX_NINPUTS = 1;

% -------------------------------------------------------------------------
% Apply the given initial condition, typically this is the only
% modification that we need to do.

ics =  [0 0 0 0 0 0 0 0 0 0];  % Initial value for each state

% -------------------------------------------------------------------------
% Load all files generated bt tge C++ program

tspan = readbin('tspan.bin', 1);
dt = tspan(2) - tspan(1);
nHorizon = length(tspan);

nX = readbin('NX.bin', MAX_NSTATES, 'int64');
nXX = prod(nX);

nU = readbin('NU.bin', MAX_NINPUTS, 'int64');
nUU = prod(nU);

u1_star_matrix = readbin('u0_star_matrix.bin', nXX);
u2_star_matrix = readbin('u1_star_matrix.bin', nXX);

descendent_matrix = readbin('descendent_matrix.bin', nXX, 'int64');

X1 = readbin('X0.bin', nX(1));
X2 = readbin('X1.bin', nX(2));
X3 = readbin('X2.bin', nX(3));
X4 = readbin('X3.bin', nX(4));
X5 = readbin('X4.bin', nX(5));
X6 = readbin('X5.bin', nX(6));
X7 = readbin('X6.bin', nX(7));
X8 = readbin('X7.bin', nX(8));
X9 = readbin('X8.bin', nX(9));
X10 = readbin('X9.bin', nX(10));

U1 = readbin('U0.bin', nU(1));
U2 = readbin('U1.bin', nU(2));

mins = readbin('state_lb.bin', MAX_NSTATES);
maxs = readbin('state_ub.bin', MAX_NSTATES);

% -------------------------------------------------------------------------
% Data holders

x1_star  = zeros(nHorizon,1);
x2_star  = zeros(nHorizon,1);
x3_star  = zeros(nHorizon,1);
x4_star  = zeros(nHorizon,1);
x5_star  = zeros(nHorizon,1);
x6_star  = zeros(nHorizon,1);
x7_star  = zeros(nHorizon,1);
x8_star  = zeros(nHorizon,1);
x9_star  = zeros(nHorizon,1);
x10_star = zeros(nHorizon,1);

u1_star  = zeros(nHorizon-1,1);
u2_star  = zeros(nHorizon-1,1);

% Start the tracing procedure
sub = cell(1,10);
for k = 1 : 10
    sub{k} = snap(ics(k), mins(k), maxs(k), nX(k)); % Start from the ICs
end

s = cell(1, 10);
ind  = sub2ind(nX, sub{:}); % Linear index for the ICs

for k = 1 : nHorizon
    [s{:}] = ind2sub(nX, ind);

    x1_star(k) = X1(s{1});
    x2_star(k) = X2(s{2});
    x3_star(k) = X3(s{3});
    x4_star(k) = X4(s{4});
    x5_star(k) = X5(s{5});
    x6_star(k) = X6(s{6});
    x7_star(k) = X7(s{7});
    x8_star(k) = X8(s{8});
    x9_star(k) = X9(s{9});
    x10_star(k) = X10(s{10});

    if k < nHorizon
        u1_star(k) = u1_star_matrix(ind, k);
        u2_star(k) = u2_star_matrix(ind, k);
    end

    ind = descendent_matrix(ind, k)+1; % C++ indexing is zero based
end

% -------------------------------------------------------------------------
% Plot all results (12 plots: 10 state-variable plots, 2 input-variable 
% plots)
% Ignore some of the plots that are not atually used!

% -------------------------------------------------------------------------
% State-variable plots

set(0,'defaultAxesFontSize', 14)
set(0,'defaultAxesFontName', 'Times')

figure
hold on
plot(tspan, x1_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(1) mins(1)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(1) maxs(1)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', "interpreter", "tex");
ylabel('\rm{x_1 (t)}', "interpreter", "tex");

figure
hold on
plot(tspan, x2_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(2) mins(2)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(2) maxs(2)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_2 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x3_star, 'b-', 'LineWidth', 2)
plot([tspan(3) tspan(end)], [mins(3) mins(3)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(3) maxs(3)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_3 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x4_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(4) mins(4)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(4) maxs(4)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_4 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x5_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(5) mins(5)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(5) maxs(5)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_5 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x6_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(6) mins(6)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(6) maxs(6)], '--r', 'LineWidth', 2)
xlabel('\rm{t}')
ylabel('\rm{x_6 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x7_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(7) mins(7)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(7) maxs(7)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_7 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x8_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(8) mins(8)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(8) maxs(8)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_8 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x9_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(9) mins(9)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(9) maxs(9)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_9 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan, x10_star, 'b-', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [mins(10) mins(10)], '--r', 'LineWidth', 2)
plot([tspan(1) tspan(end)], [maxs(10) maxs(10)], '--r', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{x_{10} (t)}', 'interpreter', 'tex');

% -------------------------------------------------------------------------
% Input-variable plots

figure
hold on
plot(tspan(1:end-1), u1_star, 'k-', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{u_1 (t)}', 'interpreter', 'tex');

figure
hold on
plot(tspan(1:end-1), u2_star, 'k-', 'LineWidth', 2)
xlabel('\rm{t}', 'interpreter', 'tex');
ylabel('\rm{u_2 (t)}', 'interpreter', 'tex');
