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
t_max = 700;
time = 1:1:t_max;
scenarios = 4;

% Change here the Scenarios to be implemented
effect_rep_rate = [1.6 2.16 3 4];                  % Effective Reproduction Rate
gamma = 1/18;                           % Average illness duration


% --------------- Simulation --------------
% Lower case value define fraction of population in every stage

% Initial values
beta_vec = gamma * effect_rep_rate;

i = randi([1 100],1)/10000;
r = 0;
e = i;   %random inition exposure
s = 1-sum([i r e]);

S_lt = zeros([t_max scenarios]);
E_lt = zeros([t_max scenarios]);
I_lt = zeros([t_max scenarios]);
R_lt = zeros([t_max scenarios]);

% Evolution
for z = 1:scenarios
    i = randi([1 100],1)/10000;              %random inition infection
    r = 0;
    e = i;  
    s = 1-sum([i r e]);
    S_lt(1,z) = s;
    E_lt(1,z) = e;
    I_lt(1,z) = i;
    R_lt(1,z) = r;

    beta = beta_vec(z);

    for j = 2:t_max
        % calculate the time derivatives
        dif_s = - beta * s * i;
        dif_e = - dif_s - sigma * e;
        dif_i = sigma * e - gamma * i;
        
        % update the states
        s = s + dif_s;
        e = e + dif_e;
        i = i + dif_i;
        r = 1 - sum([s e i]);
        
        %save at vector
        S_lt(j,z) = s;
        E_lt(j,z) = e;
        I_lt(j,z) = i;
        R_lt(j,z) = r;
    end
end

figure("Name","Infection Comparison","Color","w")
h1 = plot(I_lt(:,1));
hold on
h2 = plot(I_lt(:,2));
hold on
h3 = plot(I_lt(:,3));
hold on
h4 = plot(I_lt(:,4));
hold off

legend([h1 h2 h3 h4],["$\beta$ = " num2str(beta_vec(1))],...
                      ["$\beta$ = " num2str(beta_vec(2))],...
                      ["$\beta$ = " num2str(beta_vec(3))],...
                      ["$\beta$ = " num2str(beta_vec(4))],"Interpreter","latex")


