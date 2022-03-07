% AE6030 SPACE VEHICLE DESIGN
%
% Assignment 3: System Model Report
%
% Vinay Williams ~ K1811677
% 
% 2021
%
% SENSING FUCNTION
%
% Accepts:
%           Semimajor Axis                      [km]
%           Eccentricity                        [-]
%           Inclination                         [-]
%           Right ascension of ascending node   [deg]
%           True Anomaly                        [deg]   
%           Argument of Perigee                 [deg]
%           Simulation Time                     [days]
%           Step size                           [mins]
%           Drag Coefficient                    [-]
%           Drag Area                           [m^2]
%           Spacecraft Mass                     [kg]
%           Reflectivity Constant               []
%           Solar Radiation Pressure Area       [m^2]
%           Date                                [-]
%           Date format                         [-]
%           Tolerance                           [-]
%
% Returns:
%           Initial Date
%           Final Date
%           Initial Semimajor Axis
%           Initial Eccentricity
%           Initial Inclination 
%           Initial Argument of Perigee
%           Initial Right Ascension of Ascending Node
%           Initial True Anomaly
%           Final Semimajor Axis
%           Final Eccentricity
%           Final Inclination 
%           Final Argument of Perigee
%           Final Right Ascension of Ascending Node
%           Initial Position Magnitude (x y z)
%           Initial Velocity Magnitude (x y z)
%           Final Position Magnitude (x y z)
%           Final Velocity Magnitude (x y z)
%
% Reference : 

function [pertubations] = cowell(sma, e, inclination, raan, ta, ap, ...
                                 ndays, dtstep, cd, areadrag, scmass,...
                                 reflect, areasrp, date, dateformat, tetol)

global dtr rtd radius_earth earth_mu aunit earth_irr 
global j2 lgrav mgrav jdate0 gst0 
global rkcoef ccoef scoef ad76 bcoeff csrp0 ps
global mission

% initialize rkf78 integrator
rkcoef = 1;

% Read gravity model
gmfile = 'egm96.dat';

[ccoef, scoef] = readegm(gmfile);

% J2 coefficinets
j2 = -ccoef(3, 1);

if (strcmp(gmfile, 'egm96.dat'))
    earth_mu = 398600.4415;    
    radius_earth = 6378.1363;    
    earth_irr = 7.292115e-5;    
end


tsim = 86400.0 * ndays; % Simulation Period (seconds) 

lgrav = 18;     % Degree of gravity model
mgrav = 18;     % Tesserals in gravity model


dtstep = 60.0 * dtstep; % Step size for graphics (seconds)

% Orbital Elements
oev1(1) = sma/1000;
oev1(2) = e;
oev1(3) = inclination * dtr; 
oev1(4) = dtr * ap;
oev1(5) = dtr * raan;
oev1(6) = dtr* ta;

[fid, ad76] = read76;

% Ballistic coefficient
bcoeff = 1.0e-6 * areadrag * cd /scmass;
csrp0 = 0.000001 * reflect * ps * aunit^2 * areasrp / scmass;
    

% Initial State Vector
[ri, vi] = orb2eci(earth_mu, oev1);

% Initial position and velocity vectors
for i = 1:1:3
    yi(i) = ri(i);
    yi(i + 3) = vi(i);
end

% Initial Date (Julian)
jdate0 = juliandate(date, dateformat);


% compute initial greenwich sidereal time
gst0 = gast1(jdate0);


ti = -dtstep;
npts = 0;
    
% Initial graphics data

npts = npts + 1;
t = 0;
oev1 = eci2orb2 (earth_mu, gst0, earth_irr, t, ri, vi);

xdata(npts) = t / 86400.0;

for i = 1:1:11
        
        switch i            
            case {1, 2}                
                ydata(i, npts) = oev1(i);                
            case {3, 4, 5, 6}                
                ydata(i, npts) = rtd * oev1(i);                
            case 7                
                ydata(i, npts) = rtd * oev1(21);                
            case 8                
                ydata(i, npts) = oev1(22);                
            case 9                
                ydata(i, npts) = oev1(16);                
            case 10                
                ydata(i, npts) = rtd * oev1(15);                
            case 11                
                ydata(i, npts) = rtd * oev1(14);                
        end        
end
    


while(1)
    
   h = 30;   
        ti = ti + dtstep;       
        tf = ti + dtstep;
        
    % Integrate from inital to final time
    yfinal = rkf78('ceqm1', 6, ti, tf, h, tetol, yi);


        % Create visualisation data
        npts = npts + 1;

        % Current state vector
        for i = 1:1:3
            rf(i) = yfinal(i);
            vf(i) = yfinal(i + 3);           
        end

        % Compute current orbital elements
        oev2 = eci2orb2(earth_mu, gst0, earth_irr, tf, rf, vf);

        % Altitude check
        alt = oev2(1) * (1.0 - oev2(2)) - radius_earth;
        if (alt <= 90.0)            
            break;            
        end

        xdata(npts) = tf / 86400.0;
        
        for i = 1:1:11
     
            
            switch i                
                case {1, 2}                    
                    ydata(i, npts) = oev2(i);                    
                case {3, 4, 5, 6}                    
                    ydata(i, npts) = rtd * oev2(i);                    
                case 7                    
                    ydata(i, npts) = rtd * oev2(21);                    
                case 8                    
                    ydata(i, npts) = oev2(22);                    
                case 9                    
                    ydata(i, npts) = oev2(16);                    
                case 10                    
                    ydata(i, npts) = rtd * oev2(15);                    
                case 11                    
                    ydata(i, npts) = rtd * oev2(14);                    
            end           
        end
        
   
    yi = yfinal;

    if (tf >= tsim)        
        break;        
    end
    
end

% Final state

for i = 1:1:3
    rf(i) = yfinal(i);
    vf(i) = yfinal(i + 3);    
end

oev2 = eci2orb1(earth_mu, rf, vf);


pertubations.initial_date = datetime(jdate0,'convertfrom','juliandate');
pertubations.final_date = datetime(jdate0+ndays,'convertfrom','juliandate');

pertubations.initial_sma = oev1(1);
pertubations.initial_e= oev1(2);
pertubations.initial_i = oev1(3);
pertubations.initial_ap = oev1(4);
pertubations.initial_raan = oev1(5);
pertubations.initial_ta = oev1(6);

pertubations.initial_sma = oev2(1);
pertubations.initial_e= oev2(2);
pertubations.initial_i = oev2(3);
pertubations.initial_ap = oev2(4);
pertubations.initial_raan = oev2(5);
pertubations.initial_ta = oev2(6);

pertubations.inital_position_magnitude_x = ri(1);
pertubations.inital_position_magnitude_y = ri(2);
pertubations.inital_position_magnitude_z = ri(3);

pertubations.inital_velocity_magnitude_x = vi(1);
pertubations.inital_velocity_magnitude_y = vi(2);
pertubations.inital_velocity_magnitude_z = vi(3);

pertubations.final_position_magnitude_x = rf(1);
pertubations.final_position_magnitude_y = rf(2);
pertubations.final_position_magnitude_z = rf(3);

pertubations.final_velocity_magnitude_x = vf(1);
pertubations.final_velocity_magnitude_y = vf(2);
pertubations.final_velocity_magnitude_z = vf(3);

pertubations.xdata = xdata;
pertubations.ydata = ydata;
order = 5;
for oetype = 1:11

        switch oetype
            
            case 1
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                ylabel('Semimajor Axis (km)', 'FontSize', 12);
                title('Average Semimajor Axis evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_semimajor_axis_variation","png");
                close;
                
                pertubations.minimum_semimajor_axis = min(yFit);
                figure('visible','off')
                ylabel('Semimajor Axis (km)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Semimajor Axis evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/semimajor_axis_variation","png");
                close;
            case 2
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                ylabel('Eccentricity', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                title('Averaged Eccentricity Evolution', 'FontSize', 16);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average-eccentricity-variation","png");
                close;
                figure('visible','off')           
                ylabel('Eccentricity', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Eccentricity evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/eccentricity-variation","png");
                close;
            case 3
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                ylabel('Inclination [deg]', 'FontSize', 12);
                title('Averaged Inclination Evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average-inclination","png");
                close;
                
                pertubations.average_inclination = mean(ydata(oetype, :));
                figure('visible','off')
                ylabel('Inclination (deg)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Inclination evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/inclination","png");
                close;
            case 4
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged Argument of Perigee Evolution', 'FontSize', 16);
                ylabel('Argument of Perigee', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_semimajor_axis_variation","png");
                close;
               
                figure('visible','off')
                ylabel('Argument of Perigee (deg)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Argument of Perigee evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/ap-variation","png");
                close;
            case 5
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged RAAN Evolution', 'FontSize', 16);
                ylabel('Argument of Perigee', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_raan","png");
                close;
                figure('visible','off')
                ylabel('RAAN (deg)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('RAAN evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+ "/plots/orbit/raan-variation.png");
                close;
            case 6
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged True Anomaly Evolution', 'FontSize', 16);
                ylabel('Argument of Perigee', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_true-anomaly","png");
                close;
                figure('visible','off')
                ylabel('True Anomaly (deg)', 'FontSize', 12');
                plot(xdata, ydata(oetype, :));
                title('True Anomaly evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/true_anomaly-variation","png");
                close;
            case 7
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged Perigee Altitude Evolution', 'FontSize', 16);
                ylabel('Argument of Perigee', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_perigee_altitude_variation","png");
                close;
                figure('visible','off')           
                ylabel('Perigee Altitude (km)', 'FontSize', 12');
                plot(xdata, ydata(oetype, :));
                title('Perigee', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/perigee_altitude-variation","png");
                close;
            case 8
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged Apogee Altitude Evolution', 'FontSize', 16);
                ylabel('apogee', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_apogee_altitude_variation","png");
                close;
                
                figure('visible','off')
                ylabel('Apogee Altitude (km)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Apogee Altitude evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/apogee_altitude-variation","png");
                close;
            case 9
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged Altitude Evolution', 'FontSize', 16);
                ylabel('Altitude', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_altitude_variation","png");
                close;
                figure('visible','off')
                pertubations.average_altitude = mean(ydata(oetype, :));
                
                ylabel('Altitude (km)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Altitude evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/altitude-variation","png");
                close;
            case 10
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged East Longitude Evolution', 'FontSize', 16);
                ylabel('East Longnitude', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_east_longitude_variation","png");
                close;
                              
                figure('visible','off')
                ylabel('East Longitude (deg)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('East Longitude evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/east_longitude-variation","png");
                close;
            case 11
                figure('visible','off')
                coefficients = polyfit(xdata, ydata(oetype, :), order);
                xFit = linspace(min(xdata), max(xdata), 1000);
                yFit = polyval(coefficients , xFit);
                plot(xFit, yFit, 'r-', 'LineWidth', 2);
                grid on;
                title('Averaged Geodetic Lattitude Evolution', 'FontSize', 16);
                ylabel('Geodetic Lattitude', 'FontSize', 12);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/average_geodetic_lattitude_variation","png");
                close;
                figure('visible','off')
                ylabel('Geodetic Latitude (deg)', 'FontSize', 12);
                plot(xdata, ydata(oetype, :));
                title('Geodetic Latitude evolution', 'FontSize', 16);
                xlabel('Time (days)', 'FontSize', 12);
                grid;
                zoom on;
                saveas(gcf,"../"+mission+"/plots/orbit/geodetic_latitude-variation","png");
                close;
        end       
end