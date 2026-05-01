% Apogee distribution plot in feet AGL

clear; clc; close all;

% Given statistics in meters
mean_m   = 3552.437;
median_m = 3551.082;
std_m    = 93.255;
low_m    = 3386.271;
high_m   = 3738.432;
actual_m = 3470.65;

% Convert to feet and subtract 317 ft AGL offset
m2ft = 3.280839895;
agl_offset_ft = 317;

mean_ft   = mean_m   * m2ft - agl_offset_ft;
median_ft = median_m * m2ft - agl_offset_ft;
std_ft    = std_m    * m2ft;
low_ft    = low_m    * m2ft - agl_offset_ft;
high_ft   = high_m   * m2ft - agl_offset_ft;
actual_ft = actual_m * m2ft - agl_offset_ft;

% Create normal distribution
x = linspace(mean_ft - 4*std_ft, mean_ft + 4*std_ft, 1000);
pdf_vals = normpdf(x, mean_ft, std_ft);

figure('Color','w'); hold on; grid on; box on;

area(x, pdf_vals, 'FaceColor', [0.2 0.6 0.9], 'FaceAlpha', 0.35, ...
    'EdgeColor', [0.2 0.6 0.9], 'LineWidth', 1.5);

xline(mean_ft,   '--k', 'Mean', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(low_ft,    '--r', '95% PI Lower', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(high_ft,   '--r', '95% PI Upper', 'LineWidth', 1.8, ...
    'LabelVerticalAlignment','bottom');
xline(actual_ft,  '-m', 'Actual Apogee', 'LineWidth', 2.2, ...
    'LabelVerticalAlignment','bottom');

xlabel('Apogee (ft AGL)');
ylabel('Probability Density');
title('Apogee Distribution with 95% Prediction Interval (AGL)');
legend({'Normal PDF approximation'}, 'Location', 'northwest');

% Remove scientific notation / 10^n scaling
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat('%.0f');
ytickformat('%.4f');

text(actual_ft, max(pdf_vals)*0.85, sprintf(' Actual = %.2f ft AGL', actual_ft), ...
    'Color', 'm', 'FontWeight', 'bold');