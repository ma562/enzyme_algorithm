function [vmax_val, km_val] = M4_Algorithm_vmax_km_004_09(true_v0s, test_points, concentrations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Description
% The following function takes in the ideal measured v0 values, a vector
% of 45 possible pairs of measured v0 points and the S (concentrations)
% from the csv. It computes possible km and vmax values using two equations and two unknowns. This will
% yield 45 distinct pairs of km and vmax. We then mix up all the possible
% pairs of km and vmax to get 2025 possible pairs of km and vmax. Each km
% and vmax is passed to the SSE function to evaluate how well it represents
% the model and the best km and vmax parameters are returned to the parent
% calling function.
%
% Function Call
% [vmax_val, km_val] = M4_Algorithm_vmax_km_004_09(true_v0s, test_points, concentrations)
%
% Input Arguments
% true_v0s - a list of the measured or ideal v0s from the data
% test_points - a matrix holding all possible combinations of pairs of ideal v0
% points
% concentrations- the corresponding concentrations (S) [3.75, 7.5, 15, 30, 65, 125, 250, 500, 1000, 2000] (um)
%
% Output Arguments
% vmax_val - the value of the vmax which yields the least SSE (μM/s)
% km_val - the value of the km which yields the least SSE (μM)

SSE = realmax;      %%max possible value such that the first sse will for sure be lower than this

km = zeros(1, length(test_points(:, 1)));     %%preallocated vector of km length of 45
vmax = zeros(1, length(test_points(:, 1)));   %%preallocated vector of vmax length of 45

ind = 1;    %%index counter

%%compute possible km and vmax based on previously determined algebraic
%%equations
for i = 1:length(test_points(:, 1))
    s1 = test_points(i, 1);
    s2 = test_points(i, 3);
    v1 = test_points(i, 2);
    v2 = test_points(i, 4);
    %Improvement 1
    km(ind) = (v2 * s1 - v1 * s1) / (v1 - v2 * s1 / s2);   %%km value based on original equation (algebraic)
    vmax(ind) = v2 * (km(ind) + s2) / s2;       %%vmax value based on original equation (algebraic)
    ind = ind + 1;
end

%%if we perform all possible pairings we would get 90 * 90 = 8100 possible
%%pairs of vmax and km
%Improvement 2
for i = 1:length(km)    %%1:90
    for j = 1:length(vmax)    %%1:90
        new_SSE = M4_Algorithm_compute_SSE(true_v0s, km(i), vmax(j), concentrations);
        if(new_SSE < SSE)
            SSE = new_SSE;      %%update SSE value if current sse value is lower than previous
            vmax_val = vmax(j); %%update vmax 
            km_val = km(i);     %%update km
        end
    end
end