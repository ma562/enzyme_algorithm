function [vmax_val, km_val] = M4_Algorithm_gradient_descent_004_09(vmax, km, x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Description
% This function employs gradient descent to optimize Vmax and Km based on
% initial estimates, testing variations around these estimates to find the 
% combination that yields the lowest SSE.
%
% Function Call
% [vmax_val, km_val] = M4_Algorithm_gradient_descent_004_09(vmax, km, x, y)
%
% Input Arguments
% vmax - initial vmax value
% km - initial km value
% x - concentrations
% y - values
%
% Output Arguments
% vmax_val - the value of the vmax which yields the least SSE (μM/s)
% km_val - the value of the km which yields the least SSE (μM)

percentages = [-0.10, -0.05, 0, 0.05, 0.10];  % Variations around initial estimates
learning_rate_vmax = 0.001;
learning_rate_km = 0.01;
threshold = 1e-14;  % Small threshold for convergence
max_iterations = 20000;

% Initialize best variables
best_sse = Inf;
best_vmax = vmax;
best_km = km;

% Test variations
for pct_vmax = percentages
    for pct_km = percentages
        % Current starting points
        current_vmax = vmax * (1 + pct_vmax);
        current_km = km * (1 + pct_km);

        % Initialize for gradient descent
        working_vmax = current_vmax;
        working_km = current_km;
        previous_vmax = working_vmax;
        previous_km = working_km;

        % Perform gradient descent
        for iter = 1:max_iterations
            dSSE_dVmax = 0;
            dSSE_dKm = 0;
            for i = 1:length(x)
                pred_v = working_vmax * x(i) / (working_km + x(i));
                grad_v = 2 * (pred_v - y(i)) * (x(i) / (working_km + x(i)));
                grad_km = 2 * (pred_v - y(i)) * (-working_vmax * x(i) / ((working_km + x(i))^2));

                dSSE_dVmax = dSSE_dVmax + grad_v;
                dSSE_dKm = dSSE_dKm + grad_km;
            end

            working_vmax = working_vmax - learning_rate_vmax * dSSE_dVmax;
            working_km = working_km - learning_rate_km * dSSE_dKm;

            % Check for convergence
            if abs(working_vmax - previous_vmax) < threshold && abs(working_km - previous_km) < threshold
                break;
            end

            previous_vmax = working_vmax;
            previous_km = working_km;
        end

        % Calculate SSE for the current run
        current_sse = sum(((working_vmax .* x ./ (working_km + x)) - y).^2);

        % Update best parameters if current SSE is lower
        if current_sse < best_sse
            best_sse = current_sse;
            best_vmax = working_vmax;
            best_km = working_km;
        end
    end
end

% Return the best parameters found
vmax_val = best_vmax;
km_val = best_km;
end
