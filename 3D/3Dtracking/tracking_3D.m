function tracks3D = tracking_3D( coords3D )

%   EHarry Dec 2014


%% SETUP MOVIEINFO

nFrames = size(coords3D,1);
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),nFrames,1);
for iFrame = 1:nFrames
    coords = coords3D(iFrame).coords;
    ze = zeros(size(coords,1),2);
    movieInfo(iFrame).xCoord = [coords(:,1) ze(:,1)];
    movieInfo(iFrame).yCoord = [coords(:,2) ze(:,1)];
    movieInfo(iFrame).zCoord = [coords(:,3) ze(:,1)];
    movieInfo(iFrame).amp = ze;
end

%% GET TRACKING PARAMETERS

trackParam = trackingParameters_3D;
costMatrices = trackParam.costMatrices;
gapCloseParam = trackParam.gapCloseParam;
kalmanFunctions = trackParam.kalmanFunctions;

%% RUN TRACKER

tracks3D = trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,0);

end

