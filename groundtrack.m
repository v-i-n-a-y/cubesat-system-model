
function groundtrack(sma, e, i , raan , wp, ma, animate, speed, date, dateformat)
global mission

ae = 6371e3;  % semi-major axis of ellipsoid [m]
a = sma; % semimajor axis of a leo orbit [m]
e = e;        % eccentricity
i = i;            % inclination [deg]
Omega = raan;         % right ascension of ascending node [deg]
omega = wp;         % argument of perigee [deg]
M0 = ma;            % mean anomaly at t=0 [deg]

to = juliandate(date, dateformat);
tf = juliandate(date, dateformat) + 1;
Y = kep2cart(a,e,i,Omega,omega,M0,'mea');
step = speed;
span = 0:step:(tf-to)*86400;
option = odeset('RelTol',1e-12,'AbsTol',1e-12,'NormControl','on');
[~,r] = ode45(@Deriv,span,Y,option);
n = length(span);
ref = zeros(n,3);
lamda = zeros(n,1);
phi = zeros(n,1);
height = zeros(n,1);
for i=1:n
    to = to+step/86400;
    U = R_z(gmst(to-2400000.5));
    ref(i,1:3) = U*r(i,1:3)';
    [lamda(i),phi(i),height(i)] = Geodetic(ref(i,:));
end

figure('visible','off');
earth = imread("maps/sexymap.jpg");
    lv= size(earth,1);
    lh= size(earth,2);
    lats =  (1:lv)*180/lv - 90;
    lons =  (1:lh)*360/lh - 180;
    image(lons, lats, earth)
set(gca,'fontsize', 18);
title('TTSat Ground Track')
hold on
plot(lamda*(180/pi),phi*(180/pi),'.w')
saveas(gcf, "../"+mission+"/plots/orbit/ground-track", 'png');
close

%animation
if animate == true
    an = animatedline('Marker','*', 'MarkerSize',15, 'color', "w");
    for k = 1:n
        addpoints(an,lamda(k)*(180/pi),phi(k)*(180/pi));
        drawnow
        pause(0.2);
        clearpoints(an);
    end
end
end
%--------------------------------------------------------------------------
%
% gmst: Greenwich Mean Sidereal Time
%
% Input:
%  Mjd_UT1    Modified Julian Date UT1
%
% Output:
%  gmstime	   GMST in [rad]
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function gmstime = gmst(Mjd_UT1)
Secs = 86400;                       % Seconds per day
MJD_J2000 = 51544.5;
Mjd_0 = floor(Mjd_UT1);
UT1   = Secs*(Mjd_UT1-Mjd_0);       % [s]
T_0   = (Mjd_0  -MJD_J2000)/36525;
T     = (Mjd_UT1-MJD_J2000)/36525;
gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1...
        + (0.093104-6.2e-6*T)*T*T;  % [s]
gmstime = 2*pi*Frac(gmst/Secs);     % [rad], 0..2pi
end
%--------------------------------------------------------------------------
% 
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function res = Frac(x)
res = x-floor(x);
end
%--------------------------------------------------------------------------
%
% Geodetic: geodetic coordinates (Longitude [rad], latitude [rad],
%           altitude [m]) from given position vector (r [m])
%
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function [lon, lat, h] = Geodetic(r)
R_equ = 6378.137e3;
f     = 1/298.257223563;
epsRequ = eps*R_equ;   % Convergence criterion
e2      = f*(2-f);     % Square of eccentricity
X = r(1);              % Cartesian coordinates
Y = r(2);
Z = r(3);
rho2 = X*X + Y*Y;      % Square of distance from z-axis
% Check validity of input data
if (norm(r)==0)
    disp ( ' invalid input in Geodetic constructor\n' );
    lon = 0;
    lat = 0;
    h   = -R_equ;
    return
end
% Iteration
dZ = e2*Z;
while(1)
    ZdZ    = Z + dZ;
    Nh     = sqrt ( rho2 + ZdZ*ZdZ ); 
    SinPhi = ZdZ / Nh; % Sine of geodetic latitude
    N      = R_equ / sqrt(1-e2*SinPhi*SinPhi);
    dZ_new = N*e2*SinPhi;
    if (abs(dZ-dZ_new) < epsRequ)
        break
    end
    dZ = dZ_new;
end
% Longitude, latitude, altitude
lon = atan2(Y, X);
lat = atan2(ZdZ, sqrt(rho2));
h   = Nh - N;
end
%--------------------------------------------------------------------------
%
% Deriv: Computes the derivative of the state vector 
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function yp = Deriv(t,state)
x = state(1);
y = state(2);
z = state(3);
ae = 6371e3;             % [m]
GM = 3.986004415e14;          % [m3/s2]
J2 =  1.0826353865466185e-03; % earth's dyn. form factor (= -C20 unnormalized)
J3 = -2.5325205371813660e-06;
J4 = -1.6199892318212998e-06;
J5 = -2.2773849743222591e-07;
J6 =  5.4068710090273973e-07;
r = norm(state(1:3));
ax = -GM/r^3*x*(1 ...
   -J2*(3/2)*(ae/r)^2*(5*(z/r)^2-1));  ...
    +J3*(5/2)*(ae/r)^3*(3*(z/r)-7*(z/r)^3) ...
    -J4*(5/8)*(ae/r)^4*(3-42*(z/r)^2+63*(z/r)^4) ...
    -J5*(3/8)*(ae/r)^5*(35*(z/r)-210*(z/r)^3+231*(z/r)^5) ...
    +J6*(1/16)*(ae/r)^6*(35-945*(z/r)^2+3465*(z/r)^4-3003*(z/r)^6);
ay = y/x*ax;
az = -GM/r^3*z*(1 ...
   +J2*(3/2)*(ae/r)^2*(3-5*(z/r)^2));  ...
   +J3*(3/2)*(ae/r)^3*(10*(z/r)-(35/3)*(z/r)^3-(r/z)) ...
   -J4*(5/8)*(ae/r)^4*(15-70*(z/r)^2+63*(z/r)^4) ...
   -J5*(1/8)*(ae/r)^5*(315*(z/r)-945*(z/r)^3+693*(z/r)^5-15*(r/z)) ...
   +J6*(1/16)*(ae/r)^6*(315-2205*(z/r)^2+4851*(z/r)^4-3003*(z/r)^6);
yp = [state(4:6);ax;ay;az];
end
%--------------------------------------------------------------------------
% inputs:
%     a     semi-major axis [m]
%     ecc   eccentricity
%     inc   inclination [deg]
%     omasc right-ascension of the ascending node [deg]
%     omper argument of perigee [deg]
%     anom  time since perigee [s] or true/mean/eccentric anomaly [deg]
%     tag   string-tag: 'time', 'true' (default), 'mean' or 'eccentric'
% 
% output:
%     Y     state vector (position[m], velocity[m/s])
%
%--------------------------------------------------------------------------
function Y = kep2cart(a,ecc,inc,omasc,omper,anom,tag)
anom = anom(:);
if (nargin < 7)
    tag = 'true';
end
tag = lower(tag(1:3));
GM  = 3.986004418e14;
n   = sqrt(GM./a.^3);
if strcmp(tag,'tim')
    M = n.*anom;
    E = kepler(M,ecc,1e-10);
elseif strcmp(tag,'tru')
    f = anom/180*pi;
    E = atan2(sqrt(1-ecc.^2).*sin(f),ecc+cos(f));
elseif strcmp(tag,'ecc')
    E = anom/180*pi;
elseif strcmp(tag,'mea')
    E = kepler(anom/180*pi,ecc,1e-10);
else
    error('Undefined anomaly-type TAG.')
end
pos = [a.*(cos(E)-ecc), a.*sqrt(1-ecc.^2).*sin(E), zeros(size(anom))];
vel = [-sin(E), sqrt(1-ecc.^2).*cos(E), zeros(size(anom))];
vel = vel .* ( (n.*a./(1-ecc.*cos(E)))*[1 1 1] );
pos = rot(pos,-omper,3);
pos = rot(pos,-inc,1);
pos = rot(pos,-omasc,3);
vel = rot(vel,-omper,3);
vel = rot(vel,-inc,1);
vel = rot(vel,-omasc,3);
Y = [pos,vel]';
end
%--------------------------------------------------------------------------
% inputs:
%     mm  - mean anomaly [rad]	- matrix
%     ecc - eccentricity        - matrix
%     tol - tolerance level (def: 1e-10)
% 
% outputs:
%     ee  - eccentric anomaly [rad]
%     itr - number of iterations (optional)
%
%--------------------------------------------------------------------------
function [ee,itr] = kepler(mm,ecc,tol)
if (nargin<3)
    tol = 1e-10;
end
if ( isempty(tol) )
    error('TOL should be scalar')
end
if ( any(ecc(:) >= 1) || any(ecc(:) < 0) )
   error('Non-valid eccentricity')
end
maxitr = 20;
mm    = rem(mm,2*pi);
eeold = mm;
eenew = eeold + ecc.*sin(eeold);
itr   = 0;
while ( any(any(abs(eenew-eeold)>tol)) && (itr<maxitr) )
    itr   = itr + 1;
    eeold = eenew;
    fold  = (eeold-ecc.*sin(eeold)-mm);
    fpold = (1-ecc.*cos(eeold));
    eenew = eeold - fold./fpold;
end
if ( (itr == maxitr) && (any(any(abs(eenew-eeold)>tol))) )
   error('KEPLER didn''t achieve convergence within 20 iterations')
else
   ee = eenew;
end
end
function xout = rot(xin,alfa,i)
[r,c] = size(xin);
if ( (r ~= 3) && (c ~= 3) )
    error('XIN should be Nx3 or 3xN matrix')
elseif c == 3
    dotransp = 0;
else
    xin = xin';
    dotransp = 1;
end 
nx = size(xin,1);
na = length(alfa);
if ( (na == 1) && (nx > 1) )
    alfa = alfa*ones(nx,1);
elseif ( (na > 1) && (nx == 1) )
    xin = xin(ones(na,1),:);
elseif na ~= nx
    error('ALFA and XIN don''t match')
end
c = cos(alfa(:)/180*pi);
s = sin(alfa(:)/180*pi);
x = xin(:,1);
y = xin(:,2);
z = xin(:,3);
if (i == 1)
    xout = [    x    ,  c.*y+s.*z, -s.*y+c.*z];
elseif (i == 2)
    xout = [c.*x-s.*z,     y     ,  s.*x+c.*z];
elseif (i == 3)
    xout = [c.*x+s.*y, -s.*x+c.*y,     z     ];
else
    error('Choose axis I = 1, 2, or 3')
end
if (dotransp)
    xout = xout';
end
end
%--------------------------------------------------------------------------
%  Input:
%    angle       angle of rotation [rad]
%
%  Output:
%    rotmat      rotation matrix
%
% Last modified:   2018/01/27   M. Mahooti
%--------------------------------------------------------------------------
function rotmat = R_z(angle)
C = cos(angle);
S = sin(angle);
rotmat = zeros(3,3);
rotmat(1,1) =    C; rotmat(1,2) = S; rotmat(1,3) = 0;
rotmat(2,1) = -1*S; rotmat(2,2) = C; rotmat(2,3) = 0;
rotmat(3,1) =    0; rotmat(3,2) = 0; rotmat(3,3) = 1;
end
