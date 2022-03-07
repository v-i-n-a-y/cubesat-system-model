function oev = eci2orb2 (earth_mu, gst0, earth_irr, ut, r, v)

% convert eci state vector to complete
% set of classical orbital elements

% input

%  earth_mu    = Earth gravitational constant (km^2/sec^3)
%  gst0  = Greenwich sidereal time at 0 hours ut (radians)
%  earth_irr = Earth sidereal rotation rate (radians/second)
%  ut    = universal time (seconds)
%          (0 <= ut <= 86400)
%  r  = eci position vector (kilometers)
%  v  = eci velocity vector (kilometers/second)

% output

%  oev(1)  = semimajor axis (kilometers)
%  oev(2)  = orbital eccentricity (non-dimensional)
%            (0 <= eccentricity < 1)
%  oev(3)  = orbital inclination (radians)
%            (0 <= inclination <= pi)
%  oev(4)  = argument of perigee (radians)
%            (0 <= argument of perigee <= 2 pi)
%  oev(5)  = right ascension of ascending node (radians)
%            (0 <= raan <= 2 pi)
%  oev(6)  = true anomaly (radians)
%            (0 <= true anomaly <= 2 pi)
%  oev(7)  = orbital period (seconds)
%  oev(8)  = argument of latitude (radians)
%            (0 <= argument of latitude <= 2 pi)
%  oev(9)  = east longitude of ascending node (radians)
%            (0 <= east longitude <= 2 pi)
%  oev(10) = specific orbital energy (kilometer^2/second^2)
%  oev(11) = flight path angle (radians)
%            (-0.5 pi <= fpa <= 0.5 pi)
%  oev(12) = right ascension (radians)
%            (-2 pi <= right ascension <= 2 pi)
%  oev(13) = declination (radians)
%            (-0.5 pi <= declination <= 0.5 pi)
%  oev(14) = geodetic latitude of subpoint (radians)
%            (-0.5 pi <= latitude <= 0.5 pi)
%  oev(15) = east longitude of subpoint (radians)
%            (-2 pi <= latitude <= 2 pi)
%  oev(16) = geodetic altitude (kilometers)
%  oev(17) = geocentric radius of perigee (kilometers)
%  oev(18) = geocentric radius of apogee (kilometers)
%  oev(19) = perigee velocity (kilometers/second)
%  oev(20) = apogee velocity (kilometers/second)
%  oev(21) = geodetic altitude of perigee (kilometers)
%  oev(22) = geodetic altitude of apogee (kilometers)
%  oev(23) = geodetic latitude of perigee (radians)
%            (-0.5 pi <= latitude <= 0.5 pi)
%  oev(24) = geodetic latitude of apogee (radians)
%            (-0.5 pi <= latitude <= 0.5 pi)
%  oev(25) = mean motion (radians/second)
%  oev(26) = mean anomaly (radians)
%            (-2 pi <= mean anomaly <= 2 pi)
%  oev(27) = eccentric anomaly (radians)
%            (-2 pi <= eccentric anomaly <= 2 pi)
%  oev(28) = velocity magnitude (kilometers/second)
%  oev(29) = time until periapsis passage
%  oev(30) = time until ascending node crossing
%  oev(31) = nodal period
%  oev(32) = anomalistic period

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2.0 * pi;

% position and velocity magnitude

rmag = norm(r);

vmag = norm(v);

% position and velocity unit vectors

rhat = r / norm(r);

vhat = v / vmag;
  
% angular momentum vectors

hv = cross(r, v);

hhat = hv / norm(hv);

% eccentricity vector

vtmp = v / earth_mu;

ecc = cross(vtmp, hv);

ecc = ecc - rhat;

% semimajor axis

sma = 1 / (2 / rmag - vmag^2 / earth_mu);

% keplerian period

period = pi2 / sqrt(earth_mu / sma^3);
  
p = hhat(1) / (1 + hhat(3));

q = -hhat(2) / (1 + hhat(3));

const1 = 1 / (1 + p^2 + q^2);

fhat(1) = const1 * (1 - p^2 + q^2);
fhat(2) = const1 * 2 * p * q;
fhat(3) = -const1 * 2 * p;

ghat(1) = const1 * 2 * p * q;
ghat(2) = const1 * (1 + p^2 - q^2);
ghat(3) = const1 * 2 * q;

h = dot(ecc, ghat);

xk = dot(ecc, fhat);

x1 = dot(r, fhat);

y1 = dot(r, ghat);

% orbital eccentricity

eccm = sqrt(h * h + xk * xk);

% orbital inclination

inc = 2 * atan(sqrt(p * p + q * q));

% true longitude

xlambdat = atan3(y1, x1);

% check for equatorial orbit

if (inc > 0.0000000001)
   raan = atan3(p, q);
   xlan = mod(raan - (gst0 + earth_irr * ut), 2.0 * pi);
else
   raan = 0;
   xlan = mod(gst0 + earth_irr * ut, 2.0 * pi);
end

% check for circular orbit

if (eccm > 0.0000000001)
   argper = mod(atan3(h, xk) - raan, pi2);
else
   argper = 0;
end

% true anomaly

tanom = mod(xlambdat - raan - argper, pi2);
  
a = sin(tanom) * sqrt(1 - eccm * eccm);
b = eccm + cos(tanom);

% eccentric anomaly

eanom = atan3(a, b);

% mean anomaly

xmanom = eanom - eccm * sin(eanom);

% east longitude of subpoint

a = mod(gst0 + earth_irr * ut, pi2);

bx = sin(a);
cx = cos(a);
   
c1 = cx * r(1) + bx * r(2);

c2 = cx * r(2) - bx * r(1);

elong = atan3(c2, c1);

energy = 0.5 * vmag * vmag - earth_mu / rmag;
   
period = pi2 * sma * sqrt(sma / earth_mu);

rperigee = sma * (1 - eccm);

rapogee = sma * (1 + eccm);

vperigee = sqrt(2 * earth_mu * rapogee / (rperigee * (rapogee + rperigee)));

vapogee = sqrt(2 * earth_mu * rperigee / (rapogee * (rapogee + rperigee)));

decper = asin(sin(inc) * sin(argper));

[altper, xlatper] = geodet1(rperigee, decper);

decapo = asin(sin(inc) * sin(pi + argper));

[altapo, xlatapo] = geodet1(rapogee, decapo);

rdotv = dot(r, v);

fpa = asin(rdotv / rmag / vmag);

dec = asin(r(3) / rmag);

[alt, xlat] = geodet1(rmag, dec);

% mean motion

xmm = sqrt(earth_mu / sma * sma * sma);

ras = atan3(r(2), r(1));

if (sma > 0)
    
   % trajectory event times and orbital periods
   
   [ttpp, ttanc] = tevent (earth_mu, sma, eccm, argper, tanom);
   
   [tnodal, tanomal] = tperiod (sma, eccm, inc, argper);
   
else
    
   ttpp = 0;
   ttanc = 0;
   tnodal = 0;
   tanomal = 0;
   
end

oev(1) = sma;
oev(2) = eccm;
oev(3) = inc;
oev(4) = argper;
oev(5) = raan;
oev(6) = tanom;
oev(7) = period;
oev(8) = mod(tanom + argper, 2.0 * pi);
oev(9) = xlan;
oev(10) = energy;
oev(11) = fpa;
oev(12) = ras;
oev(13) = dec;
oev(14) = xlat;
oev(15) = elong;
oev(16) = alt;
oev(17) = rperigee;
oev(18) = rapogee;
oev(19) = vperigee;
oev(20) = vapogee;
oev(21) = altper;
oev(22) = altapo;
oev(23) = xlatper;
oev(24) = xlatapo;
oev(25) = xmm;
oev(26) = xmanom;
oev(27) = eanom;
oev(28) = vmag;
oev(29) = ttpp;
oev(30) = ttanc;
oev(31) = tnodal;
oev(32) = tanomal;



