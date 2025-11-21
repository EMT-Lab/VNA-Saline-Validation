%% VNA Validation of Saline Model ε′ and ε″ vs Measured Data

%% 1. Read Measured Data
filename = 'saline_test2.csv';
opts = detectImportOptions(filename);
opts.DataLines = [13, Inf]; % Start reading data from line 13
data = readtable(filename, opts);
data.Properties.VariableNames = {'Frequency', 'Er', 'Ei'}; % Rename columns

epsilon_0 = 8.8541878176e-12;
% Extract values
frequency = data.Frequency;         % in Hz
omega = 2 * pi * frequency;
epsilon_meas_real = data.Er;        % ε′
epsilon_meas_imag = data.Ei;        % ε″
conductivity_meas = epsilon_meas_imag .*omega .*epsilon_0;

%% 2. Input Temperature and Concentration
T = input('Enter temperature in Celsius: '); % Example: 23.7
C = 0.154;  % Physiological saline in mol/L 0.9% Saline Solution 

%% 3. Theoretical Cole-Cole Model
epsilon_0 = 8.8541878176e-12;

% Model Parameters
epsilon_static_water = 10^(1.94404 - (1.991e-3)*T);
tau_water = (3.745e-15)*(1 + (7e-5)*(T - 27.5)^2) * exp((2.2957e3) / (T + 273.15));

epsilon_s = epsilon_static_water * (1 - (3.742e-4)*T*C + 0.034*C^2 - 0.178*C + (1.515e-4)*T - (4.929e-6)*T^2);
tau = tau_water * (1.012 - (5.282e-3)*T*C + 0.032*C^2 - 0.01*C - (1.724e-3)*T + (3.766e-5)*T^2);
sigma_i = 0.174*T*C - 1.582*C^2 + 5.923*C;
alpha = (-6.348e-4)*C*T - (5.1e-2)*C^2 + (9e-2)*C;
epsilon_inf = 5.77 - 0.0274*T;

% Complex permittivity Cole-Cole Model
epsilon_complex = epsilon_inf + ...
    (epsilon_s - epsilon_inf) ./ (1 + (1j * omega * tau).^(1 - alpha)) + ...
    sigma_i ./ (1j * omega * epsilon_0);

%ε'-jε''
epsilon_model_real = real(epsilon_complex);  % ε'
epsilon_model_imag = -imag(epsilon_complex);  % ε''
conductivity =  -imag(epsilon_complex).*omega .*epsilon_0;


%% 3.5 Plot model ε' and ε'' (only theory, not measurement)
%OPTIONAL OPTIONAL OPTIONAL OPTIONAL OPTIONAL OPTIONAL OPTIONAL OPTIONAL 

figure;
plot(frequency/1e9, epsilon_model_real, 'r-', 'LineWidth', 2);
xlabel('Frequency (GHz)', 'FontSize', 14);
ylabel('\epsilon′ (Real Part)', 'FontSize', 14);
title(sprintf('Cole-Cole Model ε′ @ %.1f°C, %.3f mol/L', T, C), 'FontSize', 18);
grid on;

figure;
plot(frequency/1e9, conductivity, 'b-', 'LineWidth', 2);
xlabel('Frequency (GHz)', 'FontSize', 14);
ylabel('conductivity S/m', 'FontSize', 14);
title(sprintf('Conductivity Converted from Cole-Cole Model @ %.1f°C, %.3f mol/L', T, C), 'FontSize', 18);
grid on;

figure;
plot(frequency/1e9, epsilon_model_imag, 'b-', 'LineWidth', 2);
xlabel('Frequency (GHz)', 'FontSize', 14);
ylabel('\epsilon″ (Imaginary Part)', 'FontSize', 14);
title(sprintf('Cole-Cole Model ε″ @ %.1f°C, %.3f mol/L', T, C), 'FontSize', 18);
grid on;

%% 4. Compute Percent Error (per point)
error_Er = abs(epsilon_model_real - epsilon_meas_real) ./ abs(epsilon_model_real) * 100;
error_Ei = abs(epsilon_model_imag - epsilon_meas_imag) ./ abs(epsilon_model_imag) * 100;
error_sigma = abs(conductivity - conductivity_meas) ./ abs(conductivity) * 100;

%% 5. Report Average Errors
fprintf('\n--- Saline Model Validation @ T = %.1f°C, C = %.3f mol/L ---\n', T, C);
fprintf('Average %% Error in ε′ (real part):  %.2f%%\n', mean(error_Er));
fprintf('Average %% Error in ε″ (imag part):  %.2f%%\n', mean(error_Ei));
fprintf('Average %% Error in σ:  %.2f%%\n', mean(error_sigma));

%% 6. Plot Comparison
figure;

subplot(3,1,1);
plot(frequency/1e9, epsilon_meas_real, 'b--', 'DisplayName', 'Measured ε′'); hold on;
plot(frequency/1e9, epsilon_model_real, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Model ε′');
xlabel('Frequency (GHz)'); ylabel('\epsilon′');
title('Comparison of Real Part (\epsilon′)');
legend; grid on;

subplot(3,1,2);
plot(frequency/1e9, epsilon_meas_imag, 'b--', 'DisplayName', 'Measured ε″'); hold on;
plot(frequency/1e9, epsilon_model_imag, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Model ε″');
xlabel('Frequency (GHz)'); ylabel('\epsilon″');
title('Comparison of Imaginary Part (\epsilon″)');
legend; grid on;

subplot(3,1,3);
plot(frequency/1e9, conductivity_meas, 'b--', 'DisplayName', 'Measured ε″'); hold on;
plot(frequency/1e9, conductivity, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Model ε″');
xlabel('Frequency (GHz)'); ylabel('S/m');
title('Comparison of Conductivity');
legend; grid on;
