%% Parachute Sizing Script - Andrew Pertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    1976 Standard Atmosphere Calculator[0-1000 km]
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
drogue_Cd = 2.2;
main_Cd = 2.2;
payload_Cd = 2.2;

% Rocket
Rocket_W = 65; % rocket weight (lbs)
Payload_W = 12; % payload weight (lbs)

% Deployment
Launch_alt = 4600; % altitude of launch site, given from IREC 6.3.2 Testing and Verification Report
Drogue_alt = 10000; % drogue deployment altitude (ft)
Payload_alt = 1000; % payload deployment altitude (ft) NOTE: Same as main 

v_drogue = 65; % descent rate of drogue (ft/s)
v_payload = 25; % descent rate of payload (ft/s)

%% Calculations

%%% Get Densities
% Determine altitude density for drogue deployment
Drogue_alt = Launch_alt + Drogue_alt; % get actual alitude above sea level
Drogue_alt = convlength(Drogue_alt,'ft','km'); % convert altitude to km
[Z_drogue, Z_L_drogue, Z_U_drogue, T_drogue, P_drogue, rho_drogue,c_drogue, g_drogue, mu_drogue, nu_drogue, k_drogue, n_drogue, n_sum_drogue] = atmo(Drogue_alt,0.1,1);

% Determine altitude density for payload deployment
Payload_alt = Payload_alt + Launch_alt; % get actual altitude above sea level 
Payload_alt = convlength(Payload_alt,'ft','km'); % convert altitude to km
[Z_payload, Z_L_payload, Z_U_payload, T_payload, P_payload, rho_payload,c_payload, g_payload, mu_payload, nu_payload, k_payload, n_payload, n_sum_payload] = atmo(Payload_alt,0.1,1);

%%% Get Dynamic Pressures
v_drogue = convvel(v_drogue,'ft/s','m/s');
v_payload = convvel(v_payload,'ft/s','m/s');
q_drogue = 0.5*rho_drogue(end)*v_drogue^2; % dynamic pressure for drogue
q_payload = 0.5*rho_payload(end)*v_payload^2; % dynamic pressure for payload

%%% Determine final parachute diameters
Rocket_Drogue = convforce(Rocket_W,'lbf','N');
Payload_W = convforce(Payload_W, 'lbf', 'N');
D_drogue = sqrt((4*Rocket_Drogue)/(q_drogue*drogue_Cd*pi)); % determine drogue shoot diameter
D_payload = sqrt((4*Payload_W)/(q_payload*payload_Cd*pi)); % determine payload shoot diameter
D_main = sqrt((4*(Rocket_Drogue-Payload_W))/(q_payload*payload_Cd*pi)); % determine main shoot diameter

%% Print Statement
fprintf("The maximium drogue diameter is %1.2f (m) or %1.2f (ft)\n",D_drogue, convlength(D_drogue,'m','ft'));
fprintf("The maximium payload parachute diameter is %1.2f (m) or %1.2f (ft)\n",D_payload, convlength(D_payload,'m','ft'));
fprintf("The maximium main parachute diameter is %1.2f (m) or %1.2f (ft)\n",D_main, convlength(D_main,'m','ft'));
%% Remove Paths
rmpath("Complete 1976 Standard Atmosphere\")
disp('File Path Removed')