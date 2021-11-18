% ----------------------------------------------------------------------
%                            SIR Model
% ----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a matlab implementation of the SIR model as exposed by Sargent &
% Stachurski.
%
% It simulates the evolution of Covid based on the transition of states
% (Susceptible (S), Exposed (E), infected (I) and Removed (R).
%
% The evolution of the states is governed by the equations:
% diff_s(t)     =   -beta(t)*s(t)*i(t)
% diff_e(t)     =   beta(t)*s(t)8i(t) - sig*e(t)
% diff_r(t)      =   sigma*e(t) - gamma*i(t)
% Where :
%           beta(t): the transmission rate - rate in which individuals bump
%           into others and expose them to the virus.
%           sigma: infection rate - rate at which those who are exposed
%           becomes infected
%           gamma: recovery rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
% ----------------- Parameters -----------------
sigma = 1/5.2; % Average incubation period is of 5.2 days for covid 19
gamma = 1/18; % Average illness duration
pop_size = 2e8;
effect_rep_rate = 1.6;  %for now, it's fixed
beta = gamma * effect_rep_rate;
t_max = 700;
time = 1:1:t_max;

% --------------- Simulation --------------
% Lower case value define fraction of population in every stage

% Initial values
i = randi([1 100],1)/10000;
r = 0;
e = i;   %random inition exposure
s = 1-sum([i r e]);

S_lt = zeros([1 t_max]);
E_lt = zeros([1 t_max]);
I_lt = zeros([1 t_max]);
R_lt = zeros([1 t_max]);
S_lt(1) = s * pop_size;
E_lt(1) = e * pop_size;
I_lt(1) = i * pop_size;
R_lt(1) = r * pop_size;

% Evolution
for j =2:t_max
    S_lt(j) = s * pop_size;
    E_lt(j) = e * pop_size;
    I_lt(j) = i * pop_size;
    R_lt(j) = r * pop_size;

    dif_s = - beta * s * i;
    dif_e = - dif_s - sigma * e;
    dif_i = sigma * e - gamma * i;

    s = s + dif_s;
    e = e + dif_e;
    i = i + dif_i;
    r = 1 - sum([s e i]);
end
