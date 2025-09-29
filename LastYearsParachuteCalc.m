% Mark Stallone 
% KSU High Power Rocket Team
% Standard Atmosphere & Parachute Sizing Calculator
% 8/1/2024
%% References: Anderson- Intro to Flight 8th Edition
%%             Wikipedia (U.S. Standard Atmosphere)
%%             Previous testcode.m matlab script (I made a lot of mistakes!) 
%% NOTE: Round based on number of SigFigs in the Input/Constants sections (Application Specific) 

%% 
% Range of Altitudes to evaluate
AltFT= 0:100:20000; % from 0 to 20,000 in increments of 100 feet, change :'xxx': to +/- precision

% Allocate arrays for temperature, pressure, density, density ratio, and
% speed of sound
rho = zeros(size(AltFT));
p = zeros(size(AltFT));
temp = zeros(size(AltFT));
s_sound = zeros(size(AltFT));
density_ratio = zeros(size(AltFT));

% Calculate the values for each altitude increment, store each in array 
for i = 1:length(AltFT)
    [temp(i), p(i), rho(i), s_sound(i), density_ratio(i)] = satmFT(AltFT(i));
end

%% INPUTS 
%% Rocket
m = 65; %input("What is the mass of the EMPTY rocket in lbs? "*******THIS MEANS NO PROPELLANT*******);
m = m * 0.45359237; %Conversion factor if using lbs, comment out if using kg
%m = 4.5; %input("What is the mass of the EMPTY rocket in kg? "*******THIS MEANS NO PROPELLANT*******);
vfts_drogue = 65; %input target descent velocity in ft/s (conversion 1 ft/s to m/s is .3048 m/s)
vfts_main = 25; %input target descent velocity in ft/s (conversion 1 ft/s to m/s is .3048 m/s)
v_drogue = vfts_drogue * 0.3048; %input("What is the desired descent velocity in m/s? "); set value in desired units and multiply the variable by the constant
v_main = vfts_main * 0.3048; %input("What is the desired descent velocity in m/s? "); set value in desired units and multiply the variable by the constant
Cd_Drogue = 2.2; % Cd of Drogue parachute
Cd_Main = 2.2; % Cd of Main parachute

%% Payload
m_payload = 12; %input("What is the mass of the payload in lbs? ");
m_payload = m_payload * 0.45359237; %Conversion factor if using lbs, comment out if using kg
%m = 4.5; %input("What is the mass of the payload in kg? ");
vfts_payload = 20; %input target descent velocity in ft/s (conversion 1 ft/s to m/s is .3048 m/s)
v_payload = vfts_payload * 0.3048; %input("What is the desired descent velocity in m/s? "); set value in desired units and multiply the variable by the constant
Cd_payload = 2.2; % Cd of Drogue parachute

%% Initialization of arrays
Drogue_size_m=zeros(size(AltFT));
Main_size_m=zeros(size(AltFT));
Payload_size_m=zeros(size(AltFT));
Drogue_size_F=zeros(size(AltFT));
Drogue_size_in=zeros(size(AltFT));
Main_size_F=zeros(size(AltFT));
Main_size_in=zeros(size(AltFT));
Payload_size_F=zeros(size(AltFT));
Payload_size_in=zeros(size(AltFT));

%% Iterative Solver
for q=1:length(AltFT)
    % Gravity 
    g = 9.79324; %m/s^2 
    rho(q) = rho(q) * 515.3788; %Conversion of rho to SI Units
    Drogue_size_m(q)= sqrt((8*m*g)./(pi*rho(q)*Cd_Drogue*v_drogue^2)); %% Diameter of drogue chute (meters)
    Main_size_m(q)= sqrt((8*m*g)./(pi*rho(q)*Cd_Main*v_main^2)); %% Diameter of main chute (meters)
    Payload_size_m(q)= sqrt((8*m_payload*g)./(pi*rho(q)*Cd_payload*v_payload^2)); %% Diameter of payload chute (meters)

    % Conversions-below
    Drogue_size_F(q)= Drogue_size_m(q).*3.28084;
    Drogue_size_in(q)=Drogue_size_F(q).*12;
    Main_size_F(q)= Main_size_m(q).*3.28084;
    Main_size_in(q)=Main_size_F(q).*12;
    Payload_size_F(q)= Payload_size_m(q).*3.28084;
    Payload_size_in(q)=Payload_size_F(q).*12;
end

%% Output for ATM Data (if needed) 
% % Create a table with the results
% results_table = table(AltFT', temp', p', rho', s_sound', density_ratio', ...
%     'VariableNames', {'Altitude_ft', 'Temperature_R', 'Pressure_lbft2', 'Density_slugft3', 'Speed_of_Sound_fts', 'Density_Ratio'});
% 
% % Write the table to an Excel file
% filename = 'atmosphere_results.xlsx';
% writetable(results_table, filename, 'Sheet', 1, 'Range', 'A1');

%% Output to Command Window
% Prompt user for deployment altitudes in feet
drogue_deploy_alt = input('Enter the drogue deployment altitude in feet: ');
main_deploy_alt = input('Enter the main deployment altitude in feet: ');
payload_deploy_alt = input('Enter the payload deployment altitude in feet: ');

% Find the indices corresponding to the deployment altitudes(Basically
% figures out where the selected altitude is located in the array) 
[~, drogue_index] = min(abs(AltFT - drogue_deploy_alt));
[~, main_index] = min(abs(AltFT - main_deploy_alt));
[~, payload_index] = min(abs(AltFT - payload_deploy_alt));

% Extract the corresponding parachute sizes
drogue_size_m = Drogue_size_m(drogue_index);
drogue_size_F = Drogue_size_F(drogue_index);
drogue_size_in = Drogue_size_in(drogue_index);

main_size_m = Main_size_m(main_index);
main_size_F = Main_size_F(main_index);
main_size_in = Main_size_in(main_index);

payload_size_m = Payload_size_m(payload_index);
payload_size_F = Payload_size_F(payload_index);
payload_size_in = Payload_size_in(payload_index);

% Create a table with the results
deployment_altitudes = [drogue_deploy_alt; main_deploy_alt; payload_deploy_alt];
parachute_sizes_m = [drogue_size_m; main_size_m; payload_size_m];
parachute_sizes_F = [drogue_size_F; main_size_F; payload_size_F];
parachute_sizes_in = [drogue_size_in; main_size_in; payload_size_in];

results_table = table(deployment_altitudes, parachute_sizes_m, parachute_sizes_F, parachute_sizes_in, ...
    'VariableNames', {'Deployment_Altitude_ft', 'Parachute_Size_m', 'Parachute_Size_ft', 'Parachute_Size_in'});

% Display the table
disp(results_table);








%% Begin Function Definition

function [T, P, rhof, s_sound, density_ratio] = satmFT(h_ft)
    %% Constants (Local Variables) 
    g = 32.2; % acceleration of gravity (ft/s^2)
    T0 = 518.69; % SLS temperature (degree R)
    P0 = 2116.2; % SLS pressure (lb/ft^2)
    rho0 = 0.0023769; % SLS density (slug/ft^3)
    R = 1716; % gas constant for air (ft^2/(s^2*R))
    a = -0.00356; % temperature lapse rate (R/ft) (converted from -6.5x10^-3 K/m)
    gamma = 1.4; % constant for air
   
    %% Conditionals for Calculation
    if h_ft <= 36089
        % Troposphere layer (0 - 36089 ft, from wikipedia) 
        T = T0 + a * h_ft; % temperature calculation, Anderson pg. 120) 
        P = P0 * (T / T0) ^ (-g / (a * R)); % pressure calculation, Anderson pg. 119)
        rhof = rho0 * ((T / T0) ^ (-g / (a * R) - 1)); % density calculation, Anderson pg. 120)

    elseif h_ft <= 65617
        % Stratosphere (36089 - 65617 ft, from wikipedia)
        h_trop_ft = 36089; % transition height (ft)
        T_trop = T0 + a * h_trop_ft; % temperature at transition (R)
        P_trop = P0 * (T_trop / T0) ^ (-g / (a * R)); % pressure at transition (lb/ft^2)

        % ^ Recoded because other calculations/variables are only
        % calculated if the first 'if' statement activates. 

        T = T_trop; % temperature is constant
        P = P_trop * exp(-g * (h_ft - h_trop_ft) / (R * T_trop)); % Pressure varies
        rhof = P / (R * T); % Ideal Gas law rearranged for density calculation

    else 
        fprintf("Out of scope-go recode it idiot");
    end

    s_sound = sqrt(gamma * R * T); % calculation for the speed of sound
    density_ratio = rhof/rho0; % density ratio
end

