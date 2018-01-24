function [tracks3D, coords3D, trackMatFull] = augment3DCoords(tracks3D, coords3D, trackMatFull, coords2D, aor, angle, axisAngle)
%AUGMENT3DCOORDS Augment 3D coords by searching for additional features
%   EHarry Dec 2014

%% GET TRACKS TO PROCESS
[coordMat,idxMat] = convStruct2MatIgnoreMS(tracks3D);
nTimePoints = size(idxMat,2);
islands = getConnectedTracksWithMissingTimePoints(idxMat);

%% PROCESS TRACK ISLANDS

if ~isempty(islands)
    
    %     % get 3D tracking parameters
    %     trackParam = trackingParameters_3D;
    %     costMatrices = trackParam.costMatrices;
    %     gapCloseParam = trackParam.gapCloseParam;
    %     kalmanFunctions = trackParam.kalmanFunctions;
    
    nIslands = length(islands);
    for iIsland = 1:nIslands
        tracksToProcess = islands(iIsland).graph;
        subIndexMat = idxMat(tracksToProcess,:);
        subCoordMat = coordMat(tracksToProcess,:);
        framesToProcess = find(any(subIndexMat == 0,1));
        nFramesToProcess = length(framesToProcess);
        [potentialNewCoords,potentialNewIndexes] = deal(repmat(struct('info',[]),nFramesToProcess,1));
        [movieInfo, indexMap] = constructMovieInfo(subIndexMat, subCoordMat, 3);
        for iFrameToProcess = 1:nFramesToProcess
            frame = framesToProcess(iFrameToProcess);
            notToProcess = cell2mat(trackMatFull(:,frame));
            [potentialNewCoords(iFrameToProcess).info, potentialNewIndexes(iFrameToProcess).info] = get3DProjections(coords2D(frame,:), notToProcess);
            n = size(potentialNewCoords(iFrameToProcess).info,1);
            if n == 0
                continue
            end
            z = zeros(n,1);
            movieInfo(frame).xCoord = [movieInfo(frame).xCoord; potentialNewCoords(iFrameToProcess).info(:,1) z];
            movieInfo(frame).yCoord = [movieInfo(frame).yCoord; potentialNewCoords(iFrameToProcess).info(:,2) z];
            movieInfo(frame).zCoord = [movieInfo(frame).zCoord; potentialNewCoords(iFrameToProcess).info(:,3) z];
            movieInfo(frame).amp = [movieInfo(frame).amp; z z];
        end
        
        % track
        %trackInfo = trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,3,0,0);
        tracker_3D = tracker3D('coords',movieInfo,'aor',aor,'angle',angle,'axisAngle',axisAngle);
        tracker_3D.optimisedTracking;
        skipThisIsland = isempty(tracker_3D.tracks);
        if ~skipThisIsland
            [~,trackInfo] = convStruct2MatNoMS(tracker_3D.tracks_filtered);
            trackInfo = [trackInfo zeros(size(trackInfo,1),nTimePoints - size(trackInfo,2))];%#ok<AGROW>
        end
        delete(tracker_3D);
        clear tracker_3D
        if skipThisIsland
            continue
        end
        
        % go over each timepoint island of interest, for each track with a timepoint
        % missing through each island, find new tracks matching the track on
        % each side of the island
        for track = tracksToProcess
            timeIslands = bwconncomp(idxMat(track,:) == 0);
            for iTimeIsland = 1:timeIslands.NumObjects
                timePoints = timeIslands.PixelIdxList{iTimeIsland};
                goodTrack = true;
                newTrackIdx = [];
                if ~ismember(1,timePoints)
                    timePointToCheck = min(timePoints) - 1;
                    idxToCheck = find(idxMat(track,timePointToCheck) == indexMap(timePointToCheck).map);
                    newTrackIdx = find(trackInfo(:,timePointToCheck) == idxToCheck);
                end
                if ~ismember(nTimePoints,timePoints)
                    timePointToCheck = max(timePoints) + 1;
                    idxToCheck = find(idxMat(track,timePointToCheck) == indexMap(timePointToCheck).map);
                    newTrackIdx = [newTrackIdx find(trackInfo(:,timePointToCheck) == idxToCheck)];%#ok<AGROW>
                end
                if length(newTrackIdx) == 2
                    goodTrack = newTrackIdx(1) == newTrackIdx(2);
                    newTrackIdx = newTrackIdx(1);
                elseif isempty(newTrackIdx)
                    goodTrack = false;
                end
                
                if goodTrack
                    nNewTimePoints = length(timePoints);
                    idxToAdd = trackInfo(newTrackIdx,timePoints);
                    coordsToAdd = zeros(nNewTimePoints,3);
                    indexPairsToAdd = zeros(nNewTimePoints,2);
                    idx = 0;
                    while goodTrack && idx < nNewTimePoints
                        idx = idx + 1;
                        time = timePoints(idx);
                        id = idxToAdd(idx);
                        nOrigianlFeatures = length(indexMap(time).map);
                        goodTrack = id > nOrigianlFeatures;
                        if goodTrack
                            coordsToAdd(idx,:) = [movieInfo(time).xCoord(id,1) movieInfo(time).yCoord(id,1) movieInfo(time).zCoord(id,1)];
                            indexPairsToAdd(idx,:) = potentialNewIndexes(framesToProcess == time).info(id - nOrigianlFeatures,:);
                        end
                    end
                    
                    if goodTrack
                        idx = 0;
                        while idx < nNewTimePoints
                            idx = idx + 1;
                            time = timePoints(idx);
                            nCoordsCurrent = size(coords3D(time).coords,1) + 1;
                            coords3D(time).coords = [coords3D(time).coords; coordsToAdd(idx,:)];
                            idxMat(track,time) = nCoordsCurrent;
                            trackMatFull{track,time} = indexPairsToAdd(idx,:);
                        end
                    end
                end
            end
        end
    end
end

%% FINALLY: FOR EVERY UNASIGNED 2D FEATURE, PROJECT AGAINST ITS BEST MATCHING FEATURE IN THE OTHER SET

for iFrame = 1:nTimePoints
    ntp = cell2mat(trackMatFull(:,iFrame));
    coords3D(iFrame).coords = [coords3D(iFrame).coords; findBest2DProjections( coords2D(iFrame,:), {unique(ntp(:,1)), unique(ntp(:,2))} )];
end

end

