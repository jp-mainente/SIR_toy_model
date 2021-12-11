#=
-------------------------------------------------------------------------
                            SIR Model
-------------------------------------------------------------------------
This is a Julia implementation of the SIR model as exposed by Sargent &
Stachurski.

It simulates the evolution of Covid based on the transition of states
(Susceptible (S), Exposed (E), infected (I) and Removed (R).

The evolution of the states is governed by the equations:
diff_s(t)     =   -beta(t)*s(t)*i(t)
diff_e(t)     =   beta(t)*s(t)8i(t) - sig*e(t)
diff_r(t)      =   sigma*e(t) - gamma*i(t)
Where :
        beta(t): the transmission rate - rate in which individuals bump
        into others and expose them to the virus.
        sigma: infection rate - rate at which those who are exposed
        becomes infected
        gamma: recovery rate
-------------------------------------------------------------------------
=#

# - Adding required packages

using Pkg;
Pkg.instantiate()
using Plots

sig                 = 1/5.2;  # Average incubation period is of 5.2 days for covid 19
t_max               = 700;
time                = 1:1:t_max; 
effect_rep_rate     = [1.6 2.16 3 4];    
gamma               = 1/18; 

# --------------- Simulation --------------

beta                = effect_rep_rate * gamma;

s_lt                = zeros(t_max, size(effect_rep_rate)[2]);
e_lt                = zeros(t_max, size(effect_rep_rate)[2]);
i_lt                = zeros(t_max, size(effect_rep_rate)[2]);
r_lt                = zeros(t_max, size(effect_rep_rate)[2]);

i                   = repeat([rand(1:100)],size(effect_rep_rate)[2])/10000;
e                   = i;
r                   = repeat([0],size(effect_rep_rate)[2]);
s                   = repeat([1],size(effect_rep_rate)[2]) - (i + r + e);

s_lt[1,:]           .= s;
e_lt[1,:]           .= e;
i_lt[1,:]           .= i;
r_lt[1,:]           .= r;

for j in 2:t_max
    global s,e,i,r
    # calculate the time derivatives
    dif_s            = - beta' .* s .* i;
    dif_e            = - dif_s - sig .* e;
    dif_i            = sig .* e - gamma .* i;

    # update the states
     s              = s + dif_s;
     e              = e + dif_e;
     i              = i + dif_i;
     r              = repeat([1],size(effect_rep_rate)[2]) - (s + e + i);

    #save to vector
    s_lt[j,:]       = s;
    e_lt[j,:]       = e;
    i_lt[j,:]       = i;
    r_lt[j,:]       = r;
end

labels              = ["R = $r" for r in effect_rep_rate];
plot(i_lt,
    label = labels,
    xlabel = "t", 
    ylabel = "Infected",
    title = "Comparison Between Infection Rate")
plot(r_lt,
    label = labels,
    xlabel = "t", 
    ylabel = "Recovered",
    title = "Comparison Between Infection Rate",
    legend =:bottomright)