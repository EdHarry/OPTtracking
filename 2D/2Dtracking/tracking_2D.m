function tracks2D = tracking_2D( cellCoords , toTrack )
%TRACKING_2D simultanious 2D tracking of othogonal image series
% input: cellCoords : nFrames x 2 struct with field 'coords', an n x 2
% array, [{X,Z} , {Y,Z}] coordinates
%   EHarry Nov 2014

%% INITIALISE

if nargin < 2
    toTrack = 1:2;
end

%% SETUP MOVIEINFO

nFrames = size(cellCoords,1);
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),nFrames,2);
for iFrame = 1:nFrames
    for i = toTrack
        coords = cellCoords(iFrame,i).coords;
        ze = zeros(size(coords,1),2);
        movieInfo(iFrame,i).xCoord = [coords(:,1) ze(:,1)];
        movieInfo(iFrame,i).yCoord = [coords(:,2) ze(:,1)];
        movieInfo(iFrame,i).amp = ze;
    end
end

%% GET TRACKING PARAMETERS

trackParam = trackingParameters_2D;
costMatrices = trackParam.costMatrices;
gapCloseParam = trackParam.gapCloseParam;
kalmanFunctions = trackParam.kalmanFunctions;

%% RUN TRACKER

tracks2D = repmat(struct('trackInfo',[]),2,1);

for i = toTrack
    tracks2D(i).trackInfo = trackCloseGapsKalmanSparse(movieInfo(:,i),costMatrices,gapCloseParam,kalmanFunctions,2,0,0);
end

end

