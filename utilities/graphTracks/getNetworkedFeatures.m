function networkedFeaturesInfo = getNetworkedFeatures(trackMatFull)
%GETNETWORKEDFEATURES Summary of this function goes here
%   Detailed explanation goes here

[nTracks,nFrames] = size(trackMatFull);
networkedFeatures = repmat(struct('groups',[]),nFrames,1);
indexToGroupMap = zeros(nTracks,nFrames);

for iFrame = 1:nFrames
    feats = cell2mat(trackMatFull(:,iFrame));
    anyZeros = any(feats==0,2);
    feats(anyZeros,:) = 0;
    goodIdx = find(~anyZeros);
    feats_tmp = feats(goodIdx,:);
    nGoodIdx = length(goodIdx);
    nNodes = 2 * nGoodIdx;
    
    if nNodes == 0
        continue
    end
    
    graphMat = false(nNodes);
    graphMat(sub2ind([nNodes nNodes] ,1:nNodes ,[(nGoodIdx+1):nNodes 1:nGoodIdx])) = true;
    
    for i = 1:2
        for idx = unique(feats_tmp(:,i))'
            ind = find(feats_tmp(:,i) == idx);
            nInd = length(ind);
            for iInd = 1:(nInd-1)
                for jInd = (iInd+1):nInd
                    iNode = ind([iInd jInd]) + ((i-1)*nGoodIdx);
                    [graphMat(iNode(1),iNode(2)), graphMat(iNode(2),iNode(1))] = deal(true);
                end
            end
        end
    end
    graphGroups = floodFillGraph(graphMat);
    nGroups = length(graphGroups);
    groups = cell(nGroups,1);
    
    for iGroup = 1:nGroups
        idx = graphGroups(iGroup).graph;
        idx = idx(idx <= nGoodIdx);
        idx = goodIdx(idx);
        groups{iGroup} = feats(idx,:);
        feats(idx,1) = iGroup;
    end
    networkedFeatures(iFrame).groups = groups;
    indexToGroupMap(:,iFrame) = feats(:,1);
    
end

networkedFeaturesInfo.networkedFeatures = networkedFeatures;
networkedFeaturesInfo.indexToGroupMap = indexToGroupMap;

end

