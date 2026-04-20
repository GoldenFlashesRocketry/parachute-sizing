clc;
clear;
close all;

%% Read CSV
filename = 'RRC3TF2.csv';
opts = detectImportOptions(filename);

if any(strcmp(opts.VariableNames, 'Events'))
    opts = setvartype(opts, 'Events', 'string');
end

T = readtable(filename, opts);

time = T.Time;
altitude = T.Altitude;
velocity = T.Velocity;

%% Find apogee
[apogeeAlt, idxApogee] = max(altitude);
apogeeTime = time(idxApogee);

%% Use data after apogee only
timeDescent = time(idxApogee:end);
altDescent = altitude(idxApogee:end);
velDescent = velocity(idxApogee:end);

%% Descent rate as positive downward speed
descentRate = -velDescent;

%% Find first point at or below 1000 ft after apogee
idx1000_rel = find(altDescent <= 1000, 1, 'first');

if isempty(idx1000_rel)
    error('Rocket never descends to 1000 ft in this dataset.');
end

idx1000 = idxApogee + idx1000_rel - 1;
time1000 = time(idx1000);
alt1000 = altitude(idx1000);

%% Segment 1: apogee to 1000 ft
seg1_time = time(idxApogee:idx1000);
seg1_alt  = altitude(idxApogee:idx1000);
seg1_vel  = velocity(idxApogee:idx1000);

%% Segment 2: below 1000 ft to landing/end of data
seg2_time = time(idx1000:end);
seg2_alt  = altitude(idx1000:end);
seg2_vel  = velocity(idx1000:end);

%% Average descent rate from velocity data
% Only include actual downward motion
seg1_rate = -seg1_vel(seg1_vel < 0);
seg2_rate = -seg2_vel(seg2_vel < 0);

avgRate1_velocity = mean(seg1_rate, 'omitnan');
avgRate2_velocity = mean(seg2_rate, 'omitnan');

%% Average descent rate from altitude/time
% This is often more reliable than averaging noisy pointwise velocity
avgRate1_alttime = (seg1_alt(1) - seg1_alt(end)) / (seg1_time(end) - seg1_time(1));
avgRate2_alttime = (seg2_alt(1) - seg2_alt(end)) / (seg2_time(end) - seg2_time(1));

%% Display results
fprintf('Apogee: %.2f ft at %.2f s\n', apogeeAlt, apogeeTime);
fprintf('1000 ft crossing: %.2f ft at %.2f s\n\n', alt1000, time1000);

fprintf('Average descent rate from apogee to 1000 ft:\n');
fprintf('  Mean of velocity samples: %.2f ft/s\n', avgRate1_velocity);
fprintf('  Altitude/time average:    %.2f ft/s\n\n', avgRate1_alttime);

fprintf('Average descent rate below 1000 ft:\n');
fprintf('  Mean of velocity samples: %.2f ft/s\n', avgRate2_velocity);
fprintf('  Altitude/time average:    %.2f ft/s\n', avgRate2_alttime);

%% Plot descent rate
figure;
plot(timeDescent, descentRate, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Descent Rate (ft/s)');
title('Rocket Descent Rate vs Time');
hold on;
xline(apogeeTime, '--r', 'Apogee', 'LineWidth', 1.2);
xline(time1000, '--g', '1000 ft', 'LineWidth', 1.2);

%% Plot altitude vs time with markers
figure;
plot(time, altitude, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Altitude (ft)');
title('Rocket Altitude vs Time');
hold on;
plot(apogeeTime, apogeeAlt, 'ro', 'MarkerFaceColor', 'r');
plot(time1000, alt1000, 'go', 'MarkerFaceColor', 'g');
legend('Altitude', 'Apogee', '1000 ft', 'Location', 'best');