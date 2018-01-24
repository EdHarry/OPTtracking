function [tracks3D,coords3D,trackMatFull] = generate3DTracks(tracks2D,coords2D)
%GENERATE3DTRACKS Summary of this function goes here
%   EHarry 2014

%% GET FEATURE INDEXES OF 3D TRACKS
trackMatFull = generate3DtracksMat(tracks2D);
[nTracks,nFrames] = size(trackMatFull);

%% GET CONNECTED PROJECTION GROUPS AND GENERATE 3D COORDS
networkedFeaturesInfo = getNetworkedFeatures(trackMatFull);

tracks3DIdxMat = zeros(nTracks,nFrames);
coords3D = repmat(struct('coords',[]),nFrames,1);
for iFrame = 1:nFrames
    allIdxs = cell2mat(trackMatFull(:,iFrame));
    groups = networkedFeaturesInfo.networkedFeatures(iFrame).groups;
    indexToGroupMap = networkedFeaturesInfo.indexToGroupMap(:,iFrame);
    for iGroup = 1:length(groups)
        idxs = groups{iGroup};
        coords1 = coords2D(iFrame,1).coords(idxs(:,1),:);
        coords2 = coords2D(iFrame,2).coords(idxs(:,2),:);
        z = mean([coords1(:,2); coords2(:,2)]);
        coords = [coords1(:,1) coords2(:,1) z*ones(size(coords1,1),1)];
        nNewCoords = size(coords,1);
        nCoords = size(coords3D(iFrame).coords,1);
        coords3D(iFrame).coords = [coords3D(iFrame).coords; coords];
        idxToProcess = indexToGroupMap == iGroup;
        [~,coordPos] = ismember(allIdxs(idxToProcess,:),idxs,'rows');
        newCoordIdx = nCoords + (1:nNewCoords);
        newCoordIdx = newCoordIdx(coordPos);
        tracks3DIdxMat(idxToProcess,iFrame) = newCoordIdx';
    end
end

%% GENERATE TRACKS STRUCTURE
startTime = find(any(tracks3DIdxMat ~= 0),1);
endTime = find(any(tracks3DIdxMat ~= 0),1,'last');
tracks3DIdxMat = tracks3DIdxMat(:,startTime:endTime);
seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
[nTracks,nFeats] = size(tracks3DIdxMat);

tracks3D = repmat(struct('tracksFeatIndxCG',zeros(1,nFeats),'tracksCoordAmpCG',NaN(1,8*nFeats),'seqOfEvents',seqOfEvents),nTracks,1);

for iTrack = 1:nTracks
    startTime_ind = find(tracks3DIdxMat(iTrack,:) ~= 0,1);
    endTime_ind = find(tracks3DIdxMat(iTrack,:) ~= 0,1,'last');
    tracks3D(iTrack).tracksFeatIndxCG = tracks3DIdxMat(iTrack,startTime_ind:endTime_ind);
    tracks3D(iTrack).seqOfEvents(1,1) = startTime + startTime_ind - 1;
    tracks3D(iTrack).seqOfEvents(2,1) = endTime_ind - startTime_ind + tracks3D(iTrack).seqOfEvents(1,1);
    tracks3D(iTrack).tracksCoordAmpCG = NaN(1,8*(endTime_ind - startTime_ind + 1));
end

end

