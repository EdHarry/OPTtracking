function trackInfo = filter3DTracks( trackInfo )
%FILTER3DTRACKS Summary of this function goes here
%   Detailed explanation goes here

%% CONV TO MATRICIES 

[coordMat,idxMat] = convStruct2MatNoMS(trackInfo);

%% FILTER

% by length first
tracksToRemove = filterTracksDis3D(coordMat);
coordMat(tracksToRemove,:) = [];
idxMat(tracksToRemove,:) = [];

% then by displacements
toRemove = filterTracks3D(coordMat);
for i = 1:size(toRemove,1)
    tIdx = toRemove{i,2};
    toR = toRemove{i,1};
    for tp = toR'
        coordMat(tIdx, (tp-1)*8 + (1:8)) = NaN;
        idxMat(tIdx, tp) = 0;
    end
end

% then by length again
tracksToRemove = filterTracksDis3D(coordMat);
idxMat(tracksToRemove,:) = [];

%% GENERATE TRACKS STRUCTURE
startTime = find(any(idxMat ~= 0),1);
endTime = find(any(idxMat ~= 0),1,'last');
idxMat = idxMat(:,startTime:endTime);
seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
[nTracks,nFeats] = size(idxMat);

trackInfo = repmat(struct('tracksFeatIndxCG',zeros(1,nFeats),'tracksCoordAmpCG',NaN(1,8*nFeats),'seqOfEvents',seqOfEvents),nTracks,1);

for iTrack = 1:nTracks
    startTime_ind = find(idxMat(iTrack,:) ~= 0,1);
    endTime_ind = find(idxMat(iTrack,:) ~= 0,1,'last');
    trackInfo(iTrack).tracksFeatIndxCG = idxMat(iTrack,startTime_ind:endTime_ind);
    trackInfo(iTrack).seqOfEvents(1,1) = startTime + startTime_ind - 1;
    trackInfo(iTrack).seqOfEvents(2,1) = endTime_ind - startTime_ind + trackInfo(iTrack).seqOfEvents(1,1);
    trackInfo(iTrack).tracksCoordAmpCG = NaN(1,8*(endTime_ind - startTime_ind + 1));
end

end

