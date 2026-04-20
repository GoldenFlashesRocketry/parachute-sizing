clc;
clear;
close all;

%% Read CSV
filename = 'RRC3TF2.csv';
opts = detectImportOptions(filename);

% Make sure Events is imported as string
if any(strcmp(opts.VariableNames, 'Events'))
    opts = setvartype(opts, 'Events', 'string');
end

T = readtable(filename, opts);

time     = T.Time;
altitude = T.Altitude;
velocity = T.Velocity;
events   = T.Events;   % string array now

%% Find apogee (still useful for plotting)
[apogeeAlt, idxApogee] = max(altitude);
apogeeTime = time(idxApogee);

%% Find Drogue and Main event indices
idxDrogue = find(events == "Drogue", 1, 'first');
if isempty(idxDrogue)
    error('No "Drogue" event found in Events column.');
end

idxMain = find(events == "Main", 1, 'first');
if isempty(idxMain)
    error('No "Main" event found in Events column.');
end

% Sanity check ordering
if idxMain <= idxDrogue
    error('"Main" event occurs before or at Drogue event; check data.');
end

timeDrogue = time(idxDrogue);
altDrogue  = altitude(idxDrogue);

timeMain = time(idxMain);
altMain  = altitude(idxMain);

%% Use data after apogee for plotting descent
timeDescent = time(idxApogee:end);
altDescent  = altitude(idxApogee:end);
velDescent  = velocity(idxApogee:end);

%% Descent rate as positive downward speed
descentRate = -velDescent;

%% Segment 1: Drogue to Main
seg1_time = time(idxDrogue:idxMain);
seg1_alt  = altitude(idxDrogue:idxMain);
seg1_vel  = velocity(idxDrogue:idxMain);

%% Segment 2: Main to landing/end of data
seg2_time = time(idxMain:end);
seg2_alt  = altitude(idxMain:end);
seg2_vel  = velocity(idxMain:end);

%% Average descent rate from velocity data (positive = downward)
seg1_rate = -seg1_vel(seg1_vel < 0);
seg2_rate = -seg2_vel(seg2_vel < 0);
avgRate1_velocity = mean(seg1_rate, 'omitnan');
avgRate2_velocity = mean(seg2_rate, 'omitnan');

%% Average descent rate from altitude/time (more robust)
avgRate1_alttime = (seg1_alt(1) - seg1_alt(end)) / (seg1_time(end) - seg1_time(1));
avgRate2_alttime = (seg2_alt(1) - seg2_alt(end)) / (seg2_time(end) - seg2_time(1));

%% Display results
fprintf('Apogee: %.2f ft at %.2f s\n', apogeeAlt, apogeeTime);
fprintf('Drogue event: %.2f ft at %.2f s (index %d)\n', altDrogue, timeDrogue, idxDrogue);
fprintf('Main event:   %.2f ft at %.2f s (index %d)\n\n', altMain, timeMain, idxMain);

fprintf('Average descent rate from DROGUE to MAIN:\n');
fprintf('  Mean of velocity samples: %.2f ft/s\n', avgRate1_velocity);
fprintf('  Altitude/time average:    %.2f ft/s\n\n', avgRate1_alttime);

fprintf('Average descent rate from MAIN to landing:\n');
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
xline(timeDrogue, '--m', 'Drogue', 'LineWidth', 1.2);
xline(timeMain,   '--g', 'Main',   'LineWidth', 1.2);

%% Plot altitude vs time with markers
figure;
plot(time, altitude, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Altitude (ft)');
title('Rocket Altitude vs Time');
hold on;
plot(apogeeTime, apogeeAlt, 'ro', 'MarkerFaceColor', 'r');
plot(timeDrogue, altDrogue, 'mo', 'MarkerFaceColor', 'm');
plot(timeMain,   altMain,   'go', 'MarkerFaceColor', 'g');
legend('Altitude', 'Apogee', 'Drogue', 'Main', 'Location', 'best');