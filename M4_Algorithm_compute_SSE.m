function sse = M4_Algorithm_compute_SSE(true_v0s, km, vmax, concentration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Description 
% The following function calculates the SSE of a model determined by 
% the given parameters km and vmax. The computed model then takes in the 
% concentrations (S) and computes corresponding v0 values. 
%
% sse = M4_Algorithm_compute_SSE(true_v0s, km, vmax, concentration)
%
% Input Arguments
% true_v0s- a list of the measured or ideal v0s from the data
% km- a test km value
% vmax- a test vmax value
% concentration- the corresponding concentrations (S) [3.75, 7.5, 15, 30, 65, 125, 250, 500, 1000, 2000] (um)
%
% Output Arguments
% sse- the sse found between the ideal v0s and v0s generated from km and
% vmax

test_v0s = (vmax .* concentration) ./ (km + concentration);
sse = 0;
for i = 1:length(concentration)
    sse = sse + (test_v0s(i) - true_v0s(i))^2;
end
