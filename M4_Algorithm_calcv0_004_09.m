function v_0 = M4_Algorithm_calcv0_004_09(data, column)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% The following function takes in data and the column of interest to
% calculate v0. It takes the slope at the initial 2.5% of the data and
% returns that as v0.
%
% Function Call
% M4_Algorithm_calcv0_004_09(data, column)
%
% Input Arguments
% data- data collected from the csv
% column- column of interest in the csv 
%
% Output Arguments
% v_0- the calculated v0 values based on the initial slope of the curve (μM/s)
%
% Assignment Information
%   Assignment:     M4
%   Team ID:        004-09
%   Academic Integrity:
%     [] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers we worked with: Name, login@purdue [repeat for each]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
product = data(3:end, column);
product = product(~isnan(product)); %%notice that not all products have the same length

time = data(3:length(product) + 2, 1);
%time_interval = round(0.05 * length(time));
%Improvement 3
time_interval = round(0.025 * length(time));

t_interval = time(1: time_interval);     %%first 5 % of the time interval
p_interval = product(1: time_interval);    %%first 5 % of the product interval

p_max = p_interval(length(p_interval)); %%max product
t_max = t_interval(length(t_interval)); %%max time
p_min = p_interval(1);  %min product
t_min = t_interval(1);  %min time

v_0 = (p_max - p_min) / (t_max - t_min); %%compute slope