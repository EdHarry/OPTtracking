%% SETUP

%%%%%%
nFrames = 10;
nSpotsPerFrame = 20;
virtRad = 0.5;
initCoord = 12;
initVel = 0.1;
velVar = 0.05;
angVar = 0.05;
aor = [-10 -10];
angle = 0.9 * pi / 180;
%%%%%

%% SIMULATION

[cellCoords, originalCoords] = generate2DCoordsFrom3D('nFrames',nFrames,'nFeat',nSpotsPerFrame,'initCoord',initCoord,'initVel',initVel,'mergeProb',0,'splitProb',0,'virtRad',virtRad,'angle',angle,'aor',aor,'velVar',velVar,'angVar',angVar);

%% TRACKING

tracker = OPTtracker('coords',cellCoords,'aor',aor,'angle',angle);