function [TLE] = fetchtle(ID)

if (ID == 25544)
    TLE = { ...
        'ISS (ZARYA)'; ...
        '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'; ...
        '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'};
else
    error('ISS Placeholder only - modify code for getting TLE of other spacecraft');
end

end




