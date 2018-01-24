function [ coords3D, tracks3D ] = optimise3Dtracking( coords3D )
%OPTIMISE2DTRACKING Optimises tracking based on estimated rotations
%   EHarry Nov 2014

%% INITIALISE

% save original coords
coords_original = coords3D;

% first align by com
coords3D = alignFrames3D_com(coords3D);

% perform initial tracking
tracks3D = tracking_3D(coords3D);

% previous tracks
previousTracks = [];

%% MAIN LOOP
maxTries = 100;
count = 0;
loop = true;
%loop = false;

while count < maxTries && loop
    count = count + 1;
    
    % check previous tracks
    n = length(previousTracks);
    idx = 0;
    while loop && idx < n
        idx = idx + 1;
        if compareTracks(tracks3D, previousTracks(idx).tracks3D)
            loop = false;
        end
    end
    
    if loop
        % save trackFeature info
        previousTracks = [previousTracks; struct('tracks3D',tracks3D)];%#ok<AGROW>
        
        % align coordinates
        coords3D = alignFrames3D_rotation( coords3D, tracks3D);
        
        % re-track
        tracks3D = tracking_3D(coords3D);  
    end
end

% return original data if not converged
if loop
    disp('optimise3Dtracking; warning: tracking did not converge')
    coords3D = coords_original;
    tracks3D = [];
end

%% SUBFUNCTIONS

    function tracksEqual = compareTracks(tracks3D_1, tracks3D_2)
        [~,tracksMat_1] = convStruct2MatIgnoreMS(tracks3D_1);
        [~,tracksMat_2] = convStruct2MatIgnoreMS(tracks3D_2);
        tracksEqual = all(size(tracksMat_1) == size(tracksMat_2)) && all(ismember(tracksMat_1,tracksMat_2,'rows'));
    end

end

