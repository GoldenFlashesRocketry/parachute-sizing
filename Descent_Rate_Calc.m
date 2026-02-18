%% Parachute Sizing Script - Andrew Pertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Parachute Sizing Script
%   Author:     Andrew Pertz
%   History:    
%               
%   Input:      
%   Output:     
%   Acknowledgements:       
%   Notes:                  
%   Examples:               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
%% Add Paths
addpath("Complete 1976 Standard Atmosphere\")
disp('File Paths Added')

%% Input Parameters

% Parachutes
drogue_Cd = 1.55;
main_Cd = 2.2;
payload_Cd = 1.2;
drogue_D = 3; %ft
payload_D = 4; %ft
main_D = 8; %ft
% Rocket
Rocket_W = 50.7; % rocket weight (lbs) (from OpenRocket PR2)
Payload_W = 7.5; % payload weight (lbs)

% Deployment
Launch_alt = 2920;%4600; % altitude of launch site, given from IREC 6.3.2 Testing and Verification Report
Drogue_alt = 10000; % drogue deployment altitude (ft)
Payload_alt = 1000; % payload deployment altitude (ft) NOTE: Same as main 

delta_t = 0.1; % inflation time for parachutes
deployment_v = 23.3; % velocity at drogue deployment, taken from OpenRocket

%% Calculations 

%%% Get Densities
% Determine altitude density for drogue deployment
Drogue_alt = Launch_alt + Drogue_alt; % get actual alitude above sea level
Drogue_alt = convlength(Drogue_alt,'ft','km'); % convert altitude to km
[Z_drogue, Z_L_drogue, Z_U_drogue, T_drogue, P_drogue, rho_drogue,c_drogue, g_drogue, mu_drogue, nu_drogue, k_drogue, n_drogue, n_sum_drogue] = atmo(Drogue_alt,0.1,1);
rho_drogue = rho_drogue(end);
% Determine altitude density for payload deployment
Payload_alt = Payload_alt + Launch_alt; % get actual altitude above sea level 
Payload_alt = convlength(Payload_alt,'ft','km'); % convert altitude to km
[Z_payload, Z_L_payload, Z_U_payload, T_payload, P_payload, rho_payload,c_payload, g_payload, mu_payload, nu_payload, k_payload, n_payload, n_sum_payload] = atmo(Payload_alt,0.1,1);
rho_payload = rho_payload(end);

% Solve for Velocity
Rocket_Drogue = convforce(Rocket_W,'lbf','N');
Payload_W = convforce(Payload_W, 'lbf', 'N');
drogue_D = convlength(drogue_D,'ft','m');
main_D = convlength(main_D,'ft','m');
payload_D = convlength(payload_D,'ft','m');
V_Drogue = sqrt((8*Rocket_Drogue)/(drogue_Cd*rho_drogue*pi*drogue_D^2));
V_Main = sqrt((8*(Rocket_Drogue-Payload_W))/(main_Cd*rho_payload*pi*main_D^2));
v_payload = sqrt((8*Payload_W)/(payload_Cd*rho_payload*pi*payload_D^2));

% Snatch Forces on Drogue
a_avg_drogue = (deployment_v - V_Drogue)/delta_t;
F_avg_drogue = a_avg_drogue*(Rocket_Drogue/9.81)+Rocket_Drogue;
% Snatch Forces on Main
Rocket_Main = Rocket_Drogue - Payload_W;
a_avg_main = (V_Drogue - V_Main)/delta_t;
F_avg_main = a_avg_main*(Rocket_Main/9.81)+Rocket_Main;
% Snatch Forces on Payload
a_avg_payload = (v_payload - 0)/delta_t;
F_avg_payload = a_avg_payload*(Payload_W/9.81)+Payload_W;

%% Print Statement
disp("---Descent Rate---")
fprintf("The descent under drogue is %1.2f (m/s) or %1.2f (ft/s)\n",V_Drogue, convvel(V_Drogue,'m/s','ft/s'));
fprintf("The descent under main is %1.2f (m/s) or %1.2f (ft/s)\n",V_Main, convvel(V_Main,'m/s','ft/s'));
fprintf("The descent under payload is %1.2f (m/s) or %1.2f (ft/s)\n",v_payload, convvel(v_payload,'m/s','ft/s'));
disp("---Snatch Forces---")
fprintf("The average snatch force for Drogue is %1.2f (N)\n",F_avg_drogue);
fprintf("The average snatch force for Main is %1.2f (N)\n",F_avg_main);
fprintf("The average snatch force for Payload is %1.2f (N)\n",F_avg_payload);

%% Remove Paths
rmpath("Complete 1976 Standard Atmosphere\")
disp('File Path Removed')
