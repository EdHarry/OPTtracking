function [ cellCoords, tracks2D ] = optimise2Dtracking( cellCoords )
%OPTIMISE2DTRACKING Optimises tracking based on estimated rotations
%   EHarry Nov 2014

%% INITIALISE

% save original coords
cellCoords_original = cellCoords;

% first align by com
cellCoords = alignFrames_com(cellCoords);

% perform initial tracking
tracks2D = tracking_2D(cellCoords);

% previous tracks
previousTracks = [];

%% MAIN LOOP
maxTries = 100;
count = 0;
equalTracks = false(1,2);
loop = true;
%loop = false;

while count < maxTries && loop
    count = count + 1;
    
    % check previous tracks
    n = length(previousTracks);
    idx = 0;
    while loop && idx < n
        idx = idx + 1;
        equalTracks = compareTracks(tracks2D, previousTracks(idx).tracks2D);
        if all(equalTracks)
            loop = false;
        end
    end
    
    if loop
        % save trackFeature info
        previousTracks = [previousTracks; struct('tracks2D',tracks2D)];%#ok<AGROW>
        
        % align coordinates
        tmp = alignFrames_rotation( cellCoords, tracks2D, false, find(~equalTracks));
        cellCoords(:,~equalTracks) = tmp(:,~equalTracks);
        
        % re-track
        tmp = tracking_2D(cellCoords , find(~equalTracks));
        tracks2D(~equalTracks) = tmp(~equalTracks); 
    end
end

% return original data if not converged
if loop
    disp('optimise2Dtracking; warning: tracking did not converge')
    cellCoords = cellCoords_original;
    tracks2D = [];
end

%% SUBFUNCTIONS

    function tracksEqual = compareTracks(tracks2D_1, tracks2D_2)
        tracksEqual = false(1,2);
        [~,tracksMat_1] = convStruct2MatIgnoreMS(tracks2D_1(1).trackInfo);
        [~,tracksMat_2] = convStruct2MatIgnoreMS(tracks2D_2(1).trackInfo);
        tracksEqual(1) = all(size(tracksMat_1) == size(tracksMat_2)) && all(ismember(tracksMat_1,tracksMat_2,'rows'));
        [~,tracksMat_1] = convStruct2MatIgnoreMS(tracks2D_1(2).trackInfo);
        [~,tracksMat_2] = convStruct2MatIgnoreMS(tracks2D_2(2).trackInfo);
        tracksEqual(2) = all(size(tracksMat_1) == size(tracksMat_2)) && all(ismember(tracksMat_1,tracksMat_2,'rows'));
    end

%     function optTracking
%         tracker_2D = tracker2D('coords',cellCoords);
%         
%     end

end

