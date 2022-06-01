function M4_Main_004_09()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Description 
% The following program is the main function which loads in data from
% a csv concerning enzyme concentration values over time. The program
% calls certain functions to evaluate the data and computes the model 
% yielding the least SSE for each enzyme. It then plots the measured v0s 
% versus the predicted model to demonstrate our model's incredible
% accuracy.
%% INITIALIZATION
data = readmatrix("Data_nextGen_KEtesting_allresults.csv");
%%all_V0s = [];        %%contains 5 lists of 10 v0s where row 1 is Enzyme A row 2 is Enzyme B and so on
concentrations = data(1, 2:11);   %%extract concentrations [3.75, 7.5, 15, 30, 65, 125, 250, 500, 1000, 2000] (um)
all_V0s = zeros(5, 10);     %%preallocate dimensions for speed
v0_1 = zeros(1, 10);      %%First test
v0_2 = zeros(1, 10);      %%Test duplicate

%%PART 1: DETERMINING V0S
for i = 1:5     %%iterates through Enzyme A-E
    ind = 1;  %index of vector
    for c = (2 + (i - 1) * 20): ((2 + (i - 1)* 20) + 9)
        v0_1(ind) = M4_Algorithm_calcv0_004_09(data, c);    %%append the list of v0s to first test
        ind = ind + 1;
    end
    ind = 1;  %reset vector index
    for c = (12 + (i - 1) * 20): ((12 + (i - 1)* 20) + 9)
        v0_2(ind) = M4_Algorithm_calcv0_004_09(data, c);    %%append the list of v0s to duplicate test
        ind = ind + 1;
    end
    
    v0_avg = (v0_1 + v0_2) / 2;
    all_V0s(i, :) = v0_avg;       %%append the average of test 1 and duplicates to total v0 matrix
end

%%PART 2 DETERMINING VMAX AND KM

enzymes = ["Enzyme A", "Enzyme B", "Enzyme C", "Enzyme D", "Enzyme E"];
for row = 1:5
    v0_vals = all_V0s(row, :);          %%iterate through each list of V0s.. Enzyme A, B ...E

    point_combos = zeros(45, 4);        %%point_combos is in the format of 45 rows of (conc1, v01, conc2, v02)
    %%here we will make a vector of all possible permutations of 10 points:
    %%expecting 10 * 9 / 2 = 45
    ind = 1;
    for i = 9:-1:1
        for j = 1:i
            point_combos(ind, :) = [concentrations(10-i), v0_vals(10-i), concentrations(10 - i + j), v0_vals(10 -i + j)];
            ind = ind + 1;
        end
    end

    [vmax, km] = M4_Algorithm_vmax_km_004_09(v0_vals, point_combos, concentrations);
    sse = M4_Algorithm_compute_SSE(v0_vals, km, vmax, concentrations);
    fprintf("Enzyme %s\n", enzymes(row));
    fprintf("The vmax value is %.4f μM/s\n", vmax);
    fprintf("The km value is %.4f μM\n", km);
    fprintf("Yielding an SSE of: %.4f\n", sse);
    
    %GRAPH MEASURED V0 COMPARED TO MODEL BASED ON DETERMINED KM & VMAX
    plot(concentrations, v0_vals, "bd");
    hold on
    s = strcat(string(enzymes(row)), ' Velocity (V0s) (μM/s) vs Concentrations [P] (uM)');
    title(s);
    ylabel('Velocity (V0s) (μM/s)');
    xlabel('Concentrations [P] (uM)')
    pred_v0 = (vmax .* concentrations) ./ (km + concentrations);
    plot(concentrations, pred_v0);
    grid on
    legend("Original V0s", "Model based on calculated V-max and Km", "location", "best");
    if(row ~= 5) 
        figure;
    end
end