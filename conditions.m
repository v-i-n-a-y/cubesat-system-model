% AE6030 SPACE VEHICLE DESIGN
%
% Assignment 3: System Model Report
%
% Vinay Williams ~ K1811677
% 
% 2021
%
% Inscript Inputs
% 
% Declares variables and flag for the script 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mission = "taylor";        % Mission name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           FLAGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_input = 0;
flag_orbit_data = 1;        % Source of Orbit [1 = orbital elements (user),
                            %                  2 = tle data (use),
                            %                  3 = tle data (server)]    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inscipt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_input == 0
    
    % Mission Parameters
    
    mission_type = 2;    % Mission Type
                            %   1 - leo with no prpulsion
                            %   2 - leo with propulsion
                            %   3 - high earth
                            %   4 - planetary
                           	
                            
    % Mass Parameters
    unit_number = 6;                 % [U]  CubeSat Unit Number
    mass_payload = 1.7;              % [kg] Mass of Payload
    mass_launch = unit_number*2;     % [kg] Mass at Launch
    mass_structure = 1.7;            % [kg] Mass of Structure
    mass_margin = 0.3;               % [%]  Mass Margin
    
    % Orbital Parameters
    

    if flag_orbit_data == 1
        
        % Orbital Elements
        
        inclination = 97.3785;  % [degrees] Inclination
        altitude = 504;         % [   km  ] Altitude over main focus
        eccentricity = 0.0;     % [   -   ] Eccentricity
        ta = 0;                 % [degrees] True Anomaly
        raan = 61.225;          % [degrees] Right Ascension of Ascending Node
        ap = 0;                 % [degrees] Argument of Perigee


        % Date for orbital elements
        date = '07-Feb-2021 06:00:00';          % [   -   ] Date
        date_format = 'dd-mmm-yyyy HH:MM:SS';   % [   -   ] Format of Date
        
        % Orbit Maintainence Data
        drag_coefficient = 0.0001;  % [ - ] Drag Coefficient
        drag_area = 0.08;           % [m^2] Drag Area
        reflectivity = 0.3;         % [] Reflectivity
        areasrp = 0.018;            % [m^2] Area for solar radiation pressure
        
        % Other 
        map = "maps/map1.jpg";  % Map for groundtracks (default)
        
    elseif flag_orbit_data == 2

        TLE = { ...
        'ISS (ZARYA)'; ...
        '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'; ...
        '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'};
    
    elseif flag_orbit_data == 3 
        [tle] = fetchtledata;
    else
        disp("Error: Orbit data flag out of bounds")
        exit()
    end
        
    % Power Parameters
    
    power_payload = 12;                 % [W] Power Required by Payload 
    solar_panel_area = 0.18;            % [m^2] Solar Panel Area
    solar_panel_efficiency = 0.5;       % [   ] Solar Panel Efficiency
    solar_panel_efficiency_degradation = 0.01; % [-] 
    battery_voltage = 28;               % [V] Battery Voltage
    depth_discharge = 0.3;              % [-] Depth of Discharge for battery
    battery_specific_energy = 140;      % [Wh/kg]
  

    % Propulsion Parameters
    time = 50;              % [days] Time between orbit correction 
    h0 = altitude-15;       % [km]  Initial altitude for commissioning
    i0 = inclination-0.5;   % [deg] Initial inclination for commissioning    
    safety_factor = 1.5;    % Safety Factor to use for reserve delta v
    life_span = 5;          % [yrs] Life span of spacecraft
        

    % Thermal Parameters
    thermal_total_area = 0.16; % Total Area of Spacecraft
    thermal_albedo_area = 0.06;
    thermal_solar_area = 0.06;
    thermal_material_absorptivity = 0.3;
    thermal_material_emissivity = 0.28; 
    earth_albedo = 0.4;             % Albedo of the earth (fraction of solar flux)
    earth_infrared = 258;           % Earth's Infrared

    % Atitude Control System
    coil_area = (100)/(100^2);                   % Area of 1 coil
    acs_voltage = 12;                            % Voltage of subsystem
    centroid_dist = ((unit_number^2 *10)/2)/100; % Max distance from centroid
    
    % Communications
    frequency = 2110e6;          % [Hz] Frequency
    antenna_area = 0.004;             % [m^2] Antenna Area
    eff = 0.8;                       % Efficiency
    rx_temp = 308;                   % [K] Receiver Temperature
    bandwidth = 2000;                % [Hz] Bandwidth
    receiver_area = pi*7.3^2;
    receiver_eff = 0.75;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excel Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag_input == 1
    fprintf("ERROR: Reading from excel file is currently unsupported")
    %fprintf(" Reading inputs from excel file")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    fprintf("ERROR: Unsupport value for 'flag_input variable'")
end 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Misc

%         epoch
%         intl_designator = 25544U
%         launch_year
%         launch_number
%         launch_piece
%         catalog_number 
%         classification
%         mean_derivative
%         inclination
%         raan 
%         eccentricity
%         ap
%         mean_anomaly
%         mean_motion