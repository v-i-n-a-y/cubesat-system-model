% AE6030 SPACE VEHICLE DESIGN
%
% Assignment 3: System Model Report
%
% Vinay Williams ~ K1811677
% 
% 2021
%
% MASS BUDGET
%
% Accepts:
%           Frequency 
%           Altitude 
%           Power
%           Antenna Area
%           Efficiency of Antenna
%           Receiver Temperature
%           Bandwidth
%
% Returns: 
%           Sphere Area
%           Frequency
%           Wavelength
%           Altitude
%           Free Space Loss
%           Power
%           Efficiency
%           Antenna Area
%           Effective Antenna Area
%           Transmission Gain
%           Effective Isotropic Radiated Power
%           Received Strength
%           Receiver Gain 
%           Receiver Temperature
%           Bandwidth 
%           Noise Power 
%           Signal to noise ratio

function comms = comms (frequency, altitude, power, tx_area, tx_eff, ...
                        rx_temp, bandwidth, rx_area, rx_eff)
global c steffan_boltz

% Place frequency in struct
comms.frequency = frequency; 

% Place altitude in struct
comms.altitude = altitude;

% Place rx temp in struct
% comms.rx_temp = power/(steffan_boltz * bandwidth);
comms.rx_temp = rx_temp;

% Place power in struct
comms.power = power;

% Place tx area in struct
comms.tx_area = tx_area;

% Place tx efficiency in struct
comms.tx_efficiency = tx_eff;

% Calculate effective tx area
comms.tx_effective_area = tx_eff * tx_area;

% Place rx area in struct
comms.rx_area = rx_area;

% Place rx efficiency in struct
comms.rx_efficiency = rx_eff;

% Calculate effective rx area
comms.rx_effective_area = rx_eff * rx_area;

% Place bandwidth in struct
comms.bandwidth = bandwidth; 

% Calculate wavelngth
comms.wavelength = c/frequency;

% Calculate tx gain
comms.tx_gain = (4*pi*comms.tx_effective_area)/(comms.wavelength^2);

comms.rrx_gain = (4*pi*comms.rx_effective_area)/(comms.wavelength^2);

% Calculate effective isotropic radiated power
comms.eirp = comms.tx_gain * comms.power;

% Calculate free space loss
comms.free_space_loss = (comms.wavelength/(4*pi*comms.altitude))^2;

% Calculate received signal strength
comms.received_signal_strength = ((comms.eirp)/(4*pi*comms.altitude^2))*comms.rx_effective_area; 

% Calculate rx gain
comms.rx_gain = comms.received_signal_strength/(comms.power* comms.tx_gain* (comms.wavelength/(4*pi*comms.altitude))^2);

% Calculate noise power
comms.noise_power = steffan_boltz * comms.rx_temp * bandwidth; 

% Calculate signal to noise ratio
comms.sn = (comms.eirp/(steffan_boltz*bandwidth)*(comms.wavelength/(4*pi*comms.altitude)^2 *(comms.rx_gain/rx_temp)));

end

