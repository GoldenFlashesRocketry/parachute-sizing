% Apogee distribution plot in feet AGL
% Statistics are given in meters ASL
% Ground level is 317 meters
% Recorded apogee is given directly in feet AGL

clear; clc; close all;

% Given statistics in meters ASL
mean_m_asl   = 3552.437;
median_m_asl = 3551.082;
std_m        = 93.255;
low_m_asl    = 3386.271;
high_m_asl   = 3738.432;

% Ground elevation in meters
ground_m = 317;

% Recorded apogee in feet AGL
actual_ft_agl = 10346.63;

% Conversion
m2ft = 3.280839895013123;

% Convert statistics from ASL to AGL in meters
mean_m_agl   = mean_m_asl   - ground_m;
median_m_agl = median_m_asl - ground_m;
low_m_agl    = low_m_asl    - ground_m;
high_m_agl   = high_m_asl   - ground_m;

% Build distribution in meters AGL
x_m_agl = linspace(mean_m_agl - 4*std_m, mean_m_agl + 4*std_m, 1000);
pdf_vals = normpdf(x_m_agl, mean_m_agl, std_m);

% Convert x-axis and statistics to feet for plotting
x_ft_agl      = x_m_agl * m2ft;
mean_ft_agl   = mean_m_agl * m2ft;
median_ft_agl = median_m_agl * m2ft;
std_ft        = std_m * m2ft;
low_ft_agl    = low_m_agl * m2ft;
high_ft_agl   = high_m_agl * m2ft;

figure('Color','w'); hold on; grid on; box on;

area(x_ft_agl, pdf_vals, ...
    'FaceColor', [0.2 0.6 0.9], ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', [0.2 0.6 0.9], ...
    'LineWidth', 1.5);

xline(mean_ft_agl,   '--k', 'Mean', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(low_ft_agl,    '--r', '95% PI Lower', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(high_ft_agl,   '--r', '95% PI Upper', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(actual_ft_agl, '-m', 'Recorded Apogee', 'LineWidth', 2.2, ...
    'LabelVerticalAlignment','bottom');

xlabel('Apogee (ft AGL)');
ylabel('Probability Density');
title('Apogee Distribution with 95% Prediction Interval');
legend({'Normal PDF approximation'}, 'Location', 'northwest');

% Remove scientific notation
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat('%.0f');
ytickformat('%.4f');

text(actual_ft_agl, max(pdf_vals)*0.85, ...
    sprintf(' Recorded = %.2f ft AGL', actual_ft_agl), ...
    'Color', 'm', 'FontWeight', 'bold');

% Console output
fprintf('Mean: %.2f ft AGL\\n', mean_ft_agl);
fprintf('Median: %.2f ft AGL\\n', median_ft_agl);
fprintf('Std Dev: %.2f ft\\n', std_ft);
fprintf('95%% PI Lower: %.2f ft AGL\\n', low_ft_agl);
fprintf('95%% PI Upper: %.2f ft AGL\\n', high_ft_agl);
fprintf('Recorded Apogee: %.2f ft AGL\\n', actual_ft_agl);