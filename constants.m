% AE6030 SPACE VEHICLE DESIGN
%
% Assignment 3: System Model Report
%
% Vinay Williams ~ K1811677
% 
% 2021
%
% CONSTANTS
%
% Declares:
%             Acceleration due to gravity
%             Speed of Sound
%             Steffan Boltzmann Constant 
%             Planck Constant 
%             Speed of Light 
%             Radius of the Earth
%             Standard Gravitational Constant
%             Mass of the Earth
%             Standard Gravitational Parameter of the Earth
%             Map Path
%             Solar Flux at Earth
%             Subsystem Orbit Dry Mass Percentages
%             Weight Contingencies
%             J2 Constant
%             Mission Name
%             Radians to degrees conversion factor
%             Degress to radians conversion factor
%             Standard Gravitational Parameters of the Earth's moon
%             Standard Gravitational Parameters of the Sun
%             Earth's Inertial Rotation Rate 
%             Earth's Flatenning Factor
%             Astronomical Unit
%             Solar Constant
%             Diameter of Earth's Moon
%             Diameter of Sun
%             Power Contingency Table
%             Subsystem Power Percentage Table
%             Subsystem Power Efficiency Table
%             Solar Pressure
%             Earth's Magnetic Field Strength and its altitude

global g speed_sound steffan_boltz planck c radius_earth g_const ...
       mass_earth earth_mu map solar_flux subsystem_mass ...
       j2 mission rtd dtr moon_mu sun_mu earth_irr ...
       earth_ff aunit ps dmoon dsun  subsystem_power ...
       power_duty_cycle solar_pressure earth_magnetic;

% Earth Oblatness Gravity Coefficient
j2 = 1.08262668E-3;         

% Orbit dry mass table from ref_tables.m
subsystem_mass = []; 

% Power System Efficiencies from ref_tables.m
power_duty_cycle = [];

% Subsystem power allocations from ref_tables.m
subsystem_power = [];

% Acceleration due to Gravity   
g = 9.80665;                      

% Speed of Sound
speed_sound = 299792458;    

% Steffan Boltzmann Constant
steffan_boltz = 5.67E-08;   

% Planck's Constant 
planck = 6.626E-34;            

% Speed of Light
c = 299792458;              

% Radius of Earth
radius_earth = 6378136.3;   

% Standard Gravitational Constant
g_const = 6.67408E-11;      

% Mass of Earth
mass_earth = 5.9722E24;     

% Earth Gravitational Constant
earth_mu = 3.98618E14;      

% Map for groundtracks (default)
map = "maps/map2.png";      

% Solar Flux
solar_flux = 1360;          

% Degrees to radian factor
dtr = pi /180.0;            

% Radians to degree factor
rtd = 180.0 / pi;           

% Earth's Moon's Gravitational Constant
moon_mu = 4902.800076;      

% Sun's Gravitational Constant
sun_mu = 132712440040.944;  

% Earth's Inertial Rotation Rate
earth_irr = 7.292115486e-5; 

% Earth's Flattening Factor
earth_ff = 1.0 / 298.257;   

% Astronomical Unit
aunit = 149597870.691;      

% Mission Name (default name shown) 
mission = "Test Sat";         

% Radius of earth's moon
dmoon = 1738;         
    
% Radius of the sun
dsun = 696000;    

% Solar Constant  
ps = 0.00456;           

% Acceleration to radians
atr = dtr /3600;            

% Solar Pressure
solar_pressure = solar_flux/c; 

% Earth Magnetic Field Strength
earth_magnetic = [3.1*1e-5, 6378];
