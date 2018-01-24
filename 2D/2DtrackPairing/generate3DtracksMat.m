function trackMatFull = generate3DtracksMat(tracks2D)
%GENERATE3DTRACKSMAT Summary of this function goes here
%   EHarry Nov 2014

%% TRACK PAIRING AND CONJOINED MATRIX CONSTRUCTION

% pair tracks
trackPairs = trackPairing2D(tracks2D);
nPairs = size(trackPairs,1);

% get 2D track feature indexes and coords
idxMat = repmat(struct('mat',[]),2,1);
for i = 1:2
    [~,idxMat(i).mat] = convStruct2MatIgnoreMS(tracks2D(i).trackInfo);
end
nFrames = max(size(idxMat(1).mat,2),size(idxMat(2).mat,2));

% pad matricies
for i = 1:2
    s = size(idxMat(i).mat);
    idxMat(i).mat = [idxMat(i).mat zeros([s(1) nFrames-s(2)])];
end

% construct cell array to hold dual feature indexes
trackMatFull = cell(nPairs,nFrames);

% fill in dual array
for iPair = 1:nPairs
    idx = [idxMat(1).mat(trackPairs(iPair,1),:) ; idxMat(2).mat(trackPairs(iPair,2),:)];
    for iFrame = 1:nFrames
        trackMatFull{iPair,iFrame} = idx(:,iFrame)';
    end
end

% construct additional tracks where partner track can completely fit into a
% blank space
trackIdEnds = [size(idxMat(1).mat,1) size(idxMat(2).mat,1)];
processedPairs = [0 0];
for i = 1:2
    j = mod(i,2) + 1;
    loop = true;
    while loop
        loop = false;
        uniqueTrackIds = unique(trackPairs(:,i));
        for trackId = uniqueTrackIds'
            selectPair = find(trackPairs(:,i) == trackId);
            nSelectPair = length(selectPair);
            for iSelectPair = 1:(nSelectPair-1)
                for jSelectPair = (iSelectPair+1):nSelectPair
                    pairId = selectPair([iSelectPair jSelectPair])';
                    if ismember(pairId,processedPairs,'rows')
                        continue
                    end
                    processedPairs = [processedPairs; pairId];%#ok<AGROW>
                    
                    otherTrack = cell2mat(trackMatFull(pairId,:));
                    thisTrack = otherTrack(1,i:2:end);
                    otherTrack = otherTrack(:,j:2:end);
                    
                    % find regions of features
                    cc = bwconncomp(otherTrack(1,:) ~= 0);
                    loop2 = true;
                    iIsland = 0;
                    while iIsland < cc.NumObjects && loop2
                        iIsland = iIsland + 1;
                        if all(otherTrack(2,cc.PixelIdxList{iIsland}) == 0)
                            otherTrack(2,cc.PixelIdxList{iIsland}) = otherTrack(1,cc.PixelIdxList{iIsland});
                        else
                            loop2 = false;
                        end
                    end
                    
                    % add new track if found
                    if loop2
                        loop = true;
                        trackIdEnds(j) = trackIdEnds(j) + 1;
                        extraTrackPairs = [trackId trackIdEnds(j)];
                        if i == 2
                            extraTrackPairs = extraTrackPairs([2 1]);
                        end
                        trackPairs = [trackPairs; extraTrackPairs];%#ok<AGROW>
                        extraTrackMatFull = [thisTrack; otherTrack(2,:)];
                        if i == 2
                            extraTrackMatFull = extraTrackMatFull([2 1],:);
                        end
                        extraTrackMatFull = mat2cell(extraTrackMatFull',ones(1,nFrames),2)';
                        trackMatFull = [trackMatFull; extraTrackMatFull];%#ok<AGROW>
                    end
                end
            end
        end
    end
end
nPairs = size(trackMatFull,1);

% get unique feature index pairs, remove [0 0]
uniquePairs = unique(cell2mat(trackMatFull(:)),'rows');
uniquePairs = uniquePairs(~ismember(uniquePairs,[0 0],'rows'),:);

% make a list of indexes to unique pairs that contain a 0
uniquePairsWithZero = find(uniquePairs(:,1) == 0 | uniquePairs(:,2) == 0);

% fill in cojoined matrix with unqiue ids, make a 0 if [0 0]
trackMat = zeros(nPairs,nFrames);
for iPair = 1:nPairs
    for iFrame = 1:nFrames
        pair = trackMatFull{iPair,iFrame};
        if all(pair == [0 0])
            trackMat(iPair,iFrame) = 0;
        else
            trackMat(iPair,iFrame) = find(ismember(uniquePairs,pair,'rows'));
        end
    end
end

% remove duplicates
trackMat = unique(trackMat,'rows');

%% GENERATE EXTRA TRACKS BASED ON FEATURE ID CONFLICTS

% find conflicts
conflicts = findTrackConflicts(trackMat);

% for each conflict island, generate a sub track matrix and use it to
% generate cojoined tracks
for iConflict = 1:length(conflicts)
    confList = conflicts(iConflict).graph;
    tmpMat = trackMat(confList,:);
    extraTracks = loopCojoinedTracks(tmpMat);
    
    % add to track matrix
    trackMat = [trackMat; extraTracks];%#ok<AGROW>
end

%% REMOVE CONFLICTING TRACKS

% for each conflict island, remove the track with the most missing
% timepoints. If more than one, remove the one with the most conflicts with
% other tracks, continue until no conflics remain
loop = true;
while loop
    tmpMat = trackMat;
    tmpMat(ismember(tmpMat(:),uniquePairsWithZero)) = 0;
    
    conflicts = findTrackConflicts(tmpMat);
    toRem = [];
    
    for iConflict = 1:length(conflicts)
        confList = conflicts(iConflict).graph;
        tmpMat = trackMat(confList,:);
        nSubTracks = length(confList);
        
        nZeros = zeros(nSubTracks,1);
        for iSubTrack = 1:nSubTracks
            nZeros(iSubTrack) = sum(tmpMat(iSubTrack,:) == 0) + sum(ismember(tmpMat(iSubTrack,:),uniquePairsWithZero));
        end
        
        maxIdx = find(nZeros == max(nZeros));
        nMaxIdx = length(maxIdx);
        if nMaxIdx > 1
            
            nConfs = zeros(nMaxIdx,1);
            for iMaxIdx = 1:nMaxIdx
                n = 0;
                for iSubTrack = setdiff(1:nSubTracks,maxIdx(iMaxIdx))
                    n = n + length(findTrackConflicts(tmpMat([maxIdx(iMaxIdx) iSubTrack],:)));
                end
                nConfs(iMaxIdx) = n;
            end
            
            confMaxIdx = find(nConfs == max(nConfs),1);
            maxIdx = maxIdx(confMaxIdx);
        end
        
        toRem = [toRem; confList(maxIdx)];%#ok<AGROW>
    end
    
    if isempty(toRem)
        loop = false;
    else
        trackMat(toRem,:) = [];
    end
end

trackMat(ismember(trackMat(:),uniquePairsWithZero)) = 0;
trackMat = unique(trackMat,'rows');
trackMat(all(trackMat == 0,2),:) = [];

%% RECONSTRUCT DUAL FEATURE INDEXES

nPairs = size(trackMat,1);
trackMat = trackMat(:);
trackMatFull = zeros(nPairs*nFrames,2);
trackMatFull(trackMat ~= 0,:) = uniquePairs(trackMat(trackMat ~= 0),:);

trackMatFull = reshape(mat2cell(trackMatFull,ones(1,nPairs*nFrames),2),nPairs,nFrames);


%% SUBFUNCTIONS

    function tr = loopCojoinedTracks(inputMat)
        [tr, earlyExit_] = generateCojoinedTracks(inputMat);
        if earlyExit_
            ss = size(inputMat,1);
            if mod(ss,2) == 1
                ss = [0 (ss+1) (ss-1)] / 2;
            else
                ss = [0 ss ss] / 2;
            end
            if any(ss(2:3) < 1)
                tr = [];
                return
            end
            ss = cumsum(ss);
            tr = [];
            for iii = 2:3
                index = ss(iii-1)+1 : ss(iii);
                tr_ = loopCojoinedTracks(inputMat(index,:));
                tr = [tr; tr_];%#ok<AGROW>         
            end
        end
    end

end

