function trackIdxMat = fullIdxMatrix(track,onlySubTrackStartToConsider)
%FULLIDXMATRIX Summary of this function goes here
%   EHarry Nov 2014
% new: include all possible paths between events

if nargin < 2
    onlySubTrackStartToConsider = NaN;
end

trackStartTime = track.seqOfEvents(1,1);
[~,trackedFeatureInfo] = convStruct2MatIgnoreMS(track);
trackGraphInfo = compound2Graph(track.seqOfEvents);
startNodes = trackGraphInfo.startNodes;
endNodes = trackGraphInfo.endNodes;
splitNodes = trackGraphInfo.splitNodes;
mergeNodes = trackGraphInfo.mergeNodes;

if isnan(onlySubTrackStartToConsider)
    allNodesStart = ([startNodes; endNodes; splitNodes; mergeNodes]);
    allNodesEnd = allNodesStart;
else
    allNodesStart = (startNodes(track.seqOfEvents(startNodes,3) == onlySubTrackStartToConsider));
    allNodesEnd = ([startNodes; endNodes; splitNodes; mergeNodes]);
end

[paths, earlyExit] = findAllPaths(trackGraphInfo.trackGraph, allNodesStart, allNodesEnd);
if earlyExit
    trackIdxMat = [];
    return
end

nPaths = size(paths,1);
trackIdxMat = zeros(nPaths*4,size(trackedFeatureInfo,2));
trackIdxMat_splitExtra = [];
processedSplits = [0 0];

for iPath = 1:nPaths
    path = paths{iPath};
    for node = path
        possibleSubTracks = trackGraphInfo.subTrackIdx{node};
        time = track.seqOfEvents(node,1);
        featIdx = trackedFeatureInfo(possibleSubTracks,time);
        selectedSubTrack = unique(possibleSubTracks(featIdx ~= 0));
        while length(selectedSubTrack) > 1
            jNode = node + 1;
            if jNode > length(trackGraphInfo.subTrackIdx)
                error('fullIdxMatrix:selectedSubTrack','cannot determine "selectedSubTrack"')
            end
            if track.seqOfEvents(jNode,1) == time
                selectedSubTrack = setdiff(selectedSubTrack, trackGraphInfo.subTrackIdx{jNode});
            end
        end
        trackIdxMat((4*(iPath-1))+1,time:end) = trackedFeatureInfo(selectedSubTrack,time:end);
    end
    path = path([1 end]);
    time = track.seqOfEvents(path,1);
    isSplit = ismember(path(1),splitNodes);
    isMerge = ismember(path(2),mergeNodes);
    
    if isSplit
        splitFrom = track.seqOfEvents(path(1),4);
        time(1) = time(1) - 1;
        trackIdxMat((4*(iPath-1))+1,time(1)) = trackedFeatureInfo(splitFrom,time(1));
        if ~ismember([time(1) path(1)],processedSplits,'rows')
            trackIdxMat_splitExtra = [trackIdxMat_splitExtra; fullIdxMatrix(splitTrackHere(time(1),path(1)), splitFrom)];%#ok<AGROW>
            processedSplits = [processedSplits; time(1) path(1)];%#ok<AGROW>
        end
    end
    
    if isSplit && isMerge && diff(time) > 1
        trackIdxMat((4*(iPath-1))+2,(time(1)+1):(time(2)-1)) = trackIdxMat((4*(iPath-1))+1,(time(1)+1):(time(2)-1));
    end
    if isSplit
        trackIdxMat((4*(iPath-1))+3,(time(1)+1):time(2)) = trackIdxMat((4*(iPath-1))+1,(time(1)+1):time(2));
    end
    if isMerge
        trackIdxMat((4*(iPath-1))+4,time(1):(time(2)-1)) = trackIdxMat((4*(iPath-1))+1,time(1):(time(2)-1));
    end
end

trackIdxMat = unique([trackIdxMat; trackIdxMat_splitExtra],'rows');

allZero = all(trackIdxMat == 0,2);
trackIdxMat = trackIdxMat(~allZero,:);


    function subTrack = splitTrackHere(splitTime,splitEvent)
        subTrack = track;
        eventsToRemove = subTrack.seqOfEvents(:,1) <= splitTime & subTrack.seqOfEvents(:,2) == 2;
        subTracksToRemove = subTrack.seqOfEvents(eventsToRemove,3);
        eventsToRemove = eventsToRemove | subTrack.seqOfEvents(:,2) == 1 & ismember(subTrack.seqOfEvents(:,3),subTracksToRemove);
        eventsToConvToStarts = subTrack.seqOfEvents(:,1) <= splitTime & subTrack.seqOfEvents(:,2) == 1 & ~ismember(subTrack.seqOfEvents(:,3),subTracksToRemove);
        subTrack.seqOfEvents(eventsToConvToStarts,1) = splitTime;
        subTrack.seqOfEvents(eventsToConvToStarts,4) = NaN;
        subTrack.seqOfEvents(splitEvent,4) = NaN;
        subTrack.seqOfEvents(eventsToRemove,:) = [];
        for subTrackToRemove = sort(subTracksToRemove,1,'descend')'
            subTrack.seqOfEvents(subTrack.seqOfEvents(:,3) > subTrackToRemove,3) = subTrack.seqOfEvents(subTrack.seqOfEvents(:,3) > subTrackToRemove,3) - 1;
            subTrack.seqOfEvents(subTrack.seqOfEvents(:,4) > subTrackToRemove,4) = subTrack.seqOfEvents(subTrack.seqOfEvents(:,4) > subTrackToRemove,4) - 1;
            subTrack.tracksFeatIndxCG(subTrackToRemove,:) = [];
            subTrack.tracksCoordAmpCG(subTrackToRemove,:) = [];
        end
        timesToRemove = 1 : (splitTime - trackStartTime);
        subTrack.tracksFeatIndxCG(:,timesToRemove) = [];
        timesToRemove = 1 : (8 * (splitTime - trackStartTime));
        subTrack.tracksCoordAmpCG(:,timesToRemove) = [];
    end

end

