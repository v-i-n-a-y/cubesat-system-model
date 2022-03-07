% AE6030 SPACE VEHICLE DESIGN
%
% Assignment 3: System Model Report
%
% Vinay Williams ~ K1811677
% 
% 2021
%
% THERMALS FUCNTION
%
% Accepts:
%           Altitude 
%           Area Exposed to Solar Radiation
%           Max Distance from centroid
%           Reflectivity
%           Power of ACS system
%           Voltage of ACS system
%           Area of a coil
%
% Returns:
%           Coil Amperage
%           Maximum Force due to Solar Pressure
%           Force due to specular radiation
%           Torque due to forces
%           Magnetic Field Strength at Altitude
%           Mininum number of turns needed
%
% Reference : 
%           
%

function acs = acs(altitude, sun_area, L, q, power, voltage, coil_area)

   global solar_pressure earth_magnetic
   
   % Saving variables into structure
   acs.power = power;
   acs.voltage = voltage;
   acs.coil_area = coil_area;
   acs.sun_area = sun_area;
   acs.altitude = altitude;
   acs.centroid_dist_max = L;
   acs.reflectivity = q;
   
   % Calculate the amperage
   acs.amperage = power/voltage;
   
   % Calculate the force due to solar pressure (max condition)
   acs.force_solar_pressure_max = solar_pressure * sun_area;
   
   % Calculate the force due to specular reflection
   acs.force_specular_relfection = acs.force_solar_pressure_max * 2;
   
   % Calculate the solar torque
   acs.solar_torque = (acs.force_specular_relfection+acs. ...
                      force_solar_pressure_max) * sun_area*L*(1+q);
   
   % Calculate magentic field strength (angle removed for max condition)
   acs.magnetic_field_strength = ((earth_magnetic(1)*(earth_magnetic(2)...
                             *1000)^3)/altitude^3) *sqrt(2*sin(pi/2).^2+1);
   
   % Calculate the number of turns for one coil (max condition)
   acs.turns_min = acs.solar_torque /(acs.amperage*coil_area*acs. ...
                   magnetic_field_strength);

end

