%% VNA Validation of Saline Model ε′ and ε″ vs Measured Data (Automatic Sample Detection)
clear; clc; close all;

%% 1. Detect and Read All Measured Data Samples

[files, folderPath] = uigetfile({'*.csv;*.prn', 'VNA files'}, 'Select One or More Files', 'MultiSelect', 'on');
numSamples = numel(files);

if numSamples == 0
    error('error rename your files');
end

fprintf('Detected %d measurement samples.\n', numSamples);

all_data = cell(1, numSamples);

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

for k = 1:numSamples
    [filepath,name,ext] = fileparts(files(k));
    
    if char(ext) == '.prn'
        data = readtable(char(files(k)),"FileType","text");
    elseif char(ext) == '.csv'
        opts = detectImportOptions(char(files(k)));
        opts.DataLines = [13, Inf]; % Start reading data from line 13
        data = readtable(char(files(k)), opts);
    else
        error('File neither .prn or .csv');
    end
    
    data.Properties.VariableNames = {'Frequency', 'Er', 'Ei'}; % Rename columns
 
    all_data{k} = data;
end

   warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

% Assume all samples share the same frequency vector
frequency = all_data{1}.Frequency;
omega = 2 * pi * frequency;
epsilon_0 = 8.8541878176e-12;

% Combine and average measured ε′ and ε″ across samples
epsilon_meas_real_all = zeros(length(frequency), numSamples);
epsilon_meas_imag_all = zeros(length(frequency), numSamples);

for k = 1:numSamples
    epsilon_meas_real_all(:,k) = all_data{k}.Er;
    epsilon_meas_imag_all(:,k) = all_data{k}.Ei;
end

epsilon_meas_real = mean(epsilon_meas_real_all, 2);
epsilon_meas_imag = mean(epsilon_meas_imag_all, 2);
conductivity_meas = epsilon_meas_imag .* omega .* epsilon_0;

%% 2. Input Temperature and Concentration
T = input('Enter temperature in Celsius: '); % Example: 23.7
C = 0.154;  % Physiological saline (0.9% NaCl) in mol/L

%% 3. Theoretical Cole-Cole Model
epsilon_static_water = 10^(1.94404 - (1.991e-3)*T);
tau_water = (3.745e-15)*(1 + (7e-5)*(T - 27.5)^2) * exp((2.2957e3) / (T + 273.15));

epsilon_s = epsilon_static_water * (1 - (3.742e-4)*T*C + 0.034*C^2 - 0.178*C + ...
    (1.515e-4)*T - (4.929e-6)*T^2);
tau = tau_water * (1.012 - (5.282e-3)*T*C + 0.032*C^2 - 0.01*C - ...
    (1.724e-3)*T + (3.766e-5)*T^2);
sigma_i = 0.174*T*C - 1.582*C^2 + 5.923*C;
alpha = (-6.348e-4)*C*T - (5.1e-2)*C^2 + (9e-2)*C;
epsilon_inf = 5.77 - 0.0274*T;

epsilon_complex = epsilon_inf + ...
    (epsilon_s - epsilon_inf) ./ (1 + (1j * omega * tau).^(1 - alpha)) + ...
    sigma_i ./ (1j * omega * epsilon_0);

epsilon_model_real = real(epsilon_complex);       % ε′
epsilon_model_imag = -imag(epsilon_complex);      % ε″
conductivity = -imag(epsilon_complex) .* omega .* epsilon_0;

%% 4. Compute Percent Error (per point)
error_Er_mean = abs(epsilon_model_real - epsilon_meas_real) ./ abs(epsilon_model_real) * 100;
error_Ei_mean = abs(epsilon_model_imag - epsilon_meas_imag) ./ abs(epsilon_model_imag) * 100;
error_sigma_mean = abs(conductivity - conductivity_meas) ./ abs(conductivity) * 100;

error_Er = zeros(length(frequency),numSamples);
error_Ei = zeros(length(frequency), numSamples);

for x = 1:numSamples
    error_Er(:, x) = abs(epsilon_model_real - epsilon_meas_real_all(:,x)) ./ abs(epsilon_model_real) * 100;
    error_Ei(:, x) = abs(epsilon_model_real- epsilon_meas_real_all(:,x)) ./ abs(epsilon_model_imag) * 100;
end

%% 5. Define Upper Limit of Percent Error (per point)
limit_error_Er = zeros(1, length(frequency));
for x=1:length(frequency)
    if and(frequency(x) >= 1e9, frequency(x) <= 5e9)
        limit_error_Er(x) = 2;
    else
        limit_error_Er(x) = NaN;
    end
end

limit_error_Ei = zeros(1, length(frequency));
for x = 1:length(frequency)
    if and(frequency(x) >= 1e9, frequency(x) <= 5e9)
        limit_error_Ei(x) = 3;
    else
        limit_error_Ei(x) = NaN;
    end
end

limit_error_E = zeros(1, length(frequency));
for x = 1:length(frequency)
    if and(frequency(x) >= 1e9, frequency(x) <= 5e9)
        limit_error_E(x) = 5;
    else
        limit_error_E(x) = NaN;
    end
end

limit_error_c = zeros(1, length(frequency));
for x = 1:length(frequency)
    if and(frequency(x) >= 1e9, frequency(x) <= 5e9)
        limit_error_c(x) = 3;
    else
        limit_error_c(x) = NaN;
    end
end

%% 6. Report Average Errors
fprintf('\n--- Saline Model Validation (Average of %d Samples) @ T = %.1f°C, C = %.3f mol/L ---\n', ...
    numSamples, T, C);
fprintf('Average %% Error in ε′ (real part):  %.2f%%\n', mean(error_Er_mean));
fprintf('Average %% Error in ε″ (imag part):  %.2f%%\n', mean(error_Ei_mean));
fprintf('Average %% Error in σ:  %.2f%%\n', mean(error_sigma_mean));

%% 7. Plot Comparison (average measured vs model)

% Plot comparison
figure;
figure('Units','centimeters','Position',[2, 2, 15, 16]);

subplot(3,1,1);
plot(frequency/1e9, epsilon_model_real, 'r-', 'LineWidth', 2, 'DisplayName', 'Model ε′');
plot(frequency/1e9, epsilon_meas_real, 'b--', 'DisplayName', sprintf('Avg Measured ε′ (%d samples)', numSamples)); hold on;

xlabel('Frequency (GHz)'); ylabel('\epsilon′');
title('Comparison of Real Part (\epsilon′)');
legend; 

subplot(3,1,2);
plot(frequency/1e9, epsilon_model_imag, 'r-', 'LineWidth', 2, 'DisplayName', 'Model ε″');
plot(frequency/1e9, epsilon_meas_imag, 'b--', 'DisplayName', sprintf('Avg Measured ε″ (%d samples)', numSamples)); hold on;

xlabel('Frequency (GHz)'); ylabel('\epsilon″');
title('Comparison of Imaginary Part (\epsilon″)');
legend; 

subplot(3,1,3);
plot(frequency/1e9, conductivity, 'r-', 'LineWidth', 2, 'DisplayName', 'Model σ');
plot(frequency/1e9, conductivity_meas, 'b--', 'DisplayName', sprintf('Avg Measured σ (%d samples)', numSamples)); hold on;

xlabel('Frequency (GHz)'); ylabel('Conductivity (S/m)');
title('Comparison of Conductivity');
legend; 

% Plot calculated error
figure

subplot(4,1,1)
ymax = max(max(error_Ei_mean+error_Ei_mean,[],'all'), max(limit_error_E,[],'all'));
plot(frequency/1e9, error_Er_mean + error_Ei_mean, 'b', 'DisplayName', sprintf('Average error (%d samples)', numSamples));hold on;
plot(frequency/1e9, limit_error_E, 'r--', 'DisplayName', 'Upper limit of 5% error')
xlabel('Frequency (GHz)'); ylabel('% error');
ylim([0 ymax*1.1]);
title(sprintf('Error in ε (complex permittivity) across all %d samples', numSamples));
legend;

subplot(4,1,2)
ymax = max(max(error_Er_mean,[],'all'), max(limit_error_Er,[],'all'));
plot(frequency/1e9, error_Er_mean, 'b', 'DisplayName', sprintf('Average error (%d samples)', numSamples)); hold on;
plot(frequency/1e9, limit_error_Er, 'r--', 'DisplayName', 'Upper limit of 2% error');
xlabel('Frequency (GHz)'); ylabel('% error');
ylim([0 ymax*1.1]);
title(sprintf('Error in ε′ (real permittivity) across all %d samples', numSamples));
legend;

subplot(4,1,3)
ymax = max(max(error_Ei_mean,[],'all'), max(limit_error_Ei,[],'all'));
plot(frequency/1e9, error_Ei_mean, 'b', 'DisplayName', sprintf('Average error (%d samples)', numSamples)); hold on;
plot(frequency/1e9, limit_error_Ei, 'r--', 'DisplayName', 'Upper limit of 3% error');
xlabel('Frequency (GHz)'); ylabel('% error');
ylim([0 ymax*1.1]);
title(sprintf('Error in ε″ (imaginary permittivity) across all %d samples', numSamples));
legend;

subplot(4,1,4)
ymax = max(max(error_sigma_mean,[],'all'), max(limit_error_c,[],'all'));
plot(frequency/1e9, error_sigma_mean, 'b', 'DisplayName', sprintf('Average error (%d samples)', numSamples)); hold on;
plot(frequency/1e9, limit_error_c, 'r--', 'DisplayName', 'Upper limit of 3% error');
xlabel('Frequency (GHz)'); ylabel('% error');
ylim([0 ymax*1.1]);
title(sprintf('Error in σ (conductivity) across all %d samples', numSamples));
legend;

figure 
ymax = 2;
ymax_new = ymax;
plot(frequency/1e9, limit_error_Er, 'r--', 'DisplayName', 'Upper limit of 2% error'); hold on;
for x = 1:numSamples
    plot(frequency/1e9, error_Er(:,x), 'DisplayName', string(files(x))); hold on;
    ymax_new = max(ymax, max(error_Er(:,x)));
end
xlabel('Frequency (GHz)'); ylabel('% error');
ylim([0 ymax_new*1.1]);
title(sprintf('Per sample error in ε′ (real permittivity)', numSamples));
legend('Interpreter', 'none');

%% notes to add:
% a better way of marking the upper limit on the graph itself
% do loop to check for new files and plot
% enter temperature immediately 

