function trackInfo = splitTracks2D(trackInfo)
%SPLITTRACKS2D splits 2D tracks into all possible paths through any merges
%and splits
%   EHarry Nov 2014

[~,trackedFeatureIndx,~,numSegments] = convStruct2MatIgnoreMS(trackInfo);

nTimePoints = size(trackedFeatureIndx,2);

% first single tracks (no merging or splitting)
singleTracks = numSegments == 1;
numSegments = cumsum(numSegments);
trackFeatMat = trackedFeatureIndx(numSegments(singleTracks),:);

% now do compound tracks
compoundTracks = trackInfo(~singleTracks);
for iTrack = 1:length(compoundTracks)
    trackIdxMat_comp = fullIdxMatrix(compoundTracks(iTrack));
    s = size(trackIdxMat_comp);
    trackFeatMat = [trackFeatMat; trackIdxMat_comp zeros([s(1) nTimePoints-s(2)])];%#ok<AGROW>
end

% make tracks struct
startTime = find(any(trackFeatMat ~= 0),1);
trackFeatMat = trackFeatMat(:,startTime:end);
seqOfEvents = [startTime 1 1 NaN; nTimePoints 2 1 NaN];
[nTracks,nFeats] = size(trackFeatMat);

trackInfo = repmat(struct('tracksFeatIndxCG',zeros(1,nFeats),'tracksCoordAmpCG',NaN(1,8*nFeats),'seqOfEvents',seqOfEvents),nTracks,1);

for iTrack = 1:nTracks
    trackInfo(iTrack).tracksFeatIndxCG = trackFeatMat(iTrack,:);
end

end

