function agrav = gravity (t, y)

% first order equations of orbital motion

% N degree and M order 
% gravitational acceleration

% input

%  t = siearth_mulation time (seconds)
%  y = state vector

% output

%  agrav = ECI gravitational acceleration
%          vector (km/sec/sec)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global earth_mu radius_earth earth_irr j2 gst0

global lgrav mgrav ccoef scoef

% initialize acceleration vector

agrav = zeros(1:3);

r2 = y(1) * y(1) + y(2) * y(2) + y(3) * y(3);

r1 = sqrt(r2);

r3 = r2 * r1;

if (lgrav == 0 && mgrav == 0)
   
   % Keplerian motion

elseif (lgrav == 2 && mgrav == 0)
   
   % j2 only
   
   r5 = r2 * r3;
   
   d1 = -1.5 * j2 * radius_earth * radius_earth * earth_mu / r5;
   
   d2 = 1 - 5 * y(3) * y(3) / r2;
   
   agrav(1) = agrav(1) + y(1) * d1 * d2;
   
   agrav(2) = agrav(2) + y(2) * d1 * d2;
   
   agrav(3) = agrav(3) + y(3) * d1 * (d2 + 2);
else
   
   % user-defined degree and order gravity model

   sr2 = y(1) * y(1) + y(2) * y(2);
   
   sr1 = sqrt(sr2);
   
   sphi = y(3) / r1;
   
   phi = asin(sphi);

   % right ascension of greenwich

   pmt = mod(gst0 + earth_irr * t, 2.0 * pi);

   % east longitude of the spacecraft

   lamda = mod(atan3(y(2), y(1)) - pmt, 2.0 * pi);

   im = mgrav;

   if (mgrav < lgrav) 
      im = mgrav + 1;
   end

   p = legend(lgrav, im, sphi);

   [cn, sn, tn] = angles(mgrav, lamda, phi);

   e1 = 0;
   e2 = 0;
   e3 = 0;

   d1 = radius_earth / r1;
   d2 = 1;

   for il = 1:1:lgrav
       f1 = 0;
       f2 = 0;
       f3 = 0;
       
       il1 = il + 1;

       for im = 0:1:il
           if (im > mgrav)
              break;
           end
               
           im2 = im + 2;
           im1 = im + 1;
           
           d3 = ccoef(il1, im1) * cn(im1) + scoef(il1, im1) * sn(im1);
           
           f1 = f1 + d3 * p(il1, im1);
           
           if (im2 <= il1) 
              f2 = f2 + d3 * (p(il1, im2) - tn(im1) * p(il1, im1));
           end

           if (im2 > il1) 
              f2 = f2 - d3 * tn(im1) * p(il1, im1);
           end
           
           if (im ~= 0) 
              f3 = f3 + im * (scoef(il1, im1) * cn(im1) - ccoef(il1, im1) * sn(im1)) * p(il1, im1);
           end
       end

       d2 = d2 * d1;
       
       e1 = e1 + d2 * il1 * f1;
       
       e2 = e2 + d2 * f2;
       
       e3 = e3 + d2 * f3;
   end

   d3 = earth_mu / r1;
   
   e1 = -e1 * d3 / r1;
   
   e2 = e2 * d3;
   
   e3 = e3 * d3;
   
   d1 = e1 / r1;
   
   d2 = y(3) / (r2 * sr1) * e2;
   
   d3 = e3 / sr2;

   agrav(1) = agrav(1) + (d1 - d2) * y(1) - d3 * y(2);
   
   agrav(2) = agrav(2) + (d1 - d2) * y(2) + d3 * y(1);

   agrav(3) = agrav(3) + d1 * y(3) + sr1 * e2 / r2;

end

% complete gravity acceleration vector

for i = 1:1:3
    agrav(i) = agrav(i) - earth_mu * y(i) / r3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cn, sn, tn] = angles (m, a, b)

% earth_multiple angles

cn = zeros(21, 1);
sn = zeros(21, 1);
tn = zeros(21, 1);

cn(1) = 1;
sn(1) = 0;
tn(1) = 0;

if (m == 0)
   return;
end

cn(2) = cos(a);
sn(2) = sin(a);
tn(2) = tan(b);

if (m == 1)
   return;
end

for i = 2:1:m
    cn(i + 1) = 2 * cn(2) * cn(i) - cn(i - 1);
    sn(i + 1) = 2 * cn(2) * sn(i) - sn(i - 1);
    tn(i + 1) = tn(2) + tn(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = legend (n, m, x)

% Legendre polynomials

p(1, 1) = 1;
p(2, 1) = x;

for i = 2:1:n
    p(i + 1, 1) = ((2 * i - 1) * x * p(i, 1) -(i - 1) * p(i - 1, 1)) / i;
end

if (m == 0)
   return;
end

y = sqrt(1 - x * x);

p(2, 2) = y;

if (m == 1)
   % null
else
   for i = 2:1:m
       p(i + 1, i + 1) = (2 * i - 1) * y * p(i, i);
   end
end

for i = 2:1:n
    
    i1 = i - 1;
    
    for j = 1:1:i1
        
        if (j > m) 
           break;
        end

        p(i + 1, j + 1) = (2 * i - 1) * y * p(i, j);
            
        if (i - 2 >= j) 
           p(i + 1, j + 1) = p(i + 1, j + 1) + p(i - 1, j + 1);
        end
    end
end
