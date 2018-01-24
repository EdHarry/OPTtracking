function [extraTracks, earlyExit] = generateCojoinedTracks( idxMat )
%GENERATECOJOINEDTRACKS Summary of this function goes here
%   Detailed explanation goes here

nMaxNodes = 100;

extraTracks = idxMat;

[nTracks,nFrames] = size(idxMat);

framesToConsider = 2:(nFrames-1);

kFrame = 2;
while kFrame < nFrames
    
    confs = [];
    for jFrame = intersect(kFrame:(nFrames-1),framesToConsider)
        col = extraTracks(:,jFrame);
        col_nZ = col(col ~= 0 & ~isnan(col));
        col_unique = unique(col_nZ);
        if length(col_unique) ~= length(col_nZ)
            for idx = col_unique'
                ids = idx == col;
                if sum(ids) > 1
                    ids = find(ids);
                    nIds = length(ids);
                    %                     for i = 1:(nIds-1)
                    %                         for j = (i+1):nIds
                    %                             confs = [confs; sub2ind([nTracks nTracks],ids(i),ids(j)) jFrame];%#ok<AGROW>
                    %                         end
                    %                     end
                    for i = 1:(nIds-1)
                        confs = [confs; sub2ind([nTracks nTracks],ids(i),ids(i+1)) jFrame];%#ok<AGROW>
                    end
                end
            end
        end
    end
    
    if isempty(confs)
        break
    end
    
    confs = sortrows(confs);
    
    confsDiff = diff(confs,1,1);
    confsDiff(confsDiff(:,1) ~= 0,2) = 0;
    confs(confsDiff(:,2) == 1,:) = [];
    
    [i,j] = ind2sub([nTracks nTracks],confs(:,1));
    confs = [i j confs(:,2)];
    confs = sortrows(confs,[3 1 2]);
    
    nExtraNodes = 2 * nTracks;
    framesToConsider = unique(confs(:,3)');
    
    loop = true;
    jFrame = kFrame-1;
    previousSelect = [];
    while loop
        jFrame = jFrame + 1;
        select = ismember(confs(:,3),kFrame:jFrame);
        nNodes = size(unique([confs(select,[1 3]) ; confs(select,[2 3])],'rows'),1);
        if nNodes > nMaxNodes && any(select)
            loop = false;
            if ~isempty(previousSelect)
                select = previousSelect;
                jFrame = jFrame - 1;
                nNodes = nNodesPrevious;
            end
        elseif jFrame == nFrames-1
            loop = false;
        elseif any(select)
            previousSelect = select;
            nNodesPrevious = nNodes;
        end
    end
    
    nNodes = nNodes + nExtraNodes;
    confs = confs(select,:);
    kFrame = jFrame+1;
    
    graphMat = false(nNodes);
    trackMap = zeros(nNodes,2);
    
    previousFrames = [];
    idx = 0;
    buildTransMat = false;
    for iFrame = unique(confs(:,3))'
        confs_sub = confs(confs(:,3)==iFrame,1:2);
        
        for trans = confs_sub'
            for trackId = trans'
                trackIds = [trackId iFrame];
                if ~ismember(trackIds,trackMap,'rows')
                    idx = idx + 1;
                    trackMap(idx,:) = trackIds;
                end
            end
            thisFrameIdx = trackMap(:,2) == iFrame;
            iNode = trackMap(:,1) == trans(1) & thisFrameIdx;
            jNode = trackMap(:,1) == trans(2) & thisFrameIdx;
            graphMat(iNode,jNode) = true;
            graphMat(jNode,iNode) = true;
        end
        
        if buildTransMat
            for trackId = unique(confs_sub(:))'
                thisFrameIdx = trackMap(:,2) == iFrame;
                jNode = trackMap(:,1) == trackId & thisFrameIdx;
                if any(jNode)
                    loop = true;
                    idx2 = 0;
                    n = length(previousFrames);
                    while idx2 < n && loop
                        idx2 = idx2 + 1;
                        previousFrame = previousFrames(idx2);
                        previousFrameIdx = trackMap(:,2) == previousFrame;
                        iNode = trackMap(:,1) == trackId & previousFrameIdx;
                        if any(iNode)
                            loop = false;
                            graphMat(iNode,jNode) = true;
                        end
                    end
                end
            end
        else
            buildTransMat = true;
        end
        
        previousFrames = [previousFrames iFrame];%#ok<AGROW>
    end
    
    trackMap(idx+1:idx+nTracks,:) = [(1:nTracks)' ones(nTracks,1)];
    trackMap(idx+nTracks+1:end,:) = [(1:nTracks)' nFrames*ones(nTracks,1)];
    
    startTime = [1 nFrames];
    dir = [1 -1];
    for i = 1:2
        for iTrack = 1:nTracks
            thisTrackIdx = trackMap(:,1) == iTrack;
            iNode = trackMap(:,2) == startTime(i) & thisTrackIdx;
            linked = false;
            idx2 = 0;
            while ~linked
                idx2 = idx2 + 1;
                nextTP = startTime(i) + (dir(i) * idx2);
                jNode = trackMap(:,2) == nextTP & thisTrackIdx;
                if any(jNode)
                    if i == 2
                        tmp = iNode;
                        iNode = jNode;
                        jNode = tmp;
                    end
                    graphMat(iNode,jNode) = true;
                    linked = true;
                end
            end
        end
    end
    
    [paths, earlyExit] = findAllPaths(graphMat, ((idx+1):(idx+nTracks))', ((idx+nTracks+1):nNodes)');
    if earlyExit
        extraTracks = [];
        return
    end
    
    nPaths = size(paths,1);
    
    extraTracks_ = zeros(nPaths,nFrames);
    
    for iPath = 1:nPaths
        path = paths{iPath};
        for node = path
            map = trackMap(node,:);
            extraTracks_(iPath,map(2):end) = extraTracks(map(1),map(2):end);
        end
    end
    
    extraTracks = unique(extraTracks_,'rows');
    nTracks = size(extraTracks,1);
    
end

extraTracks = extraTracks(~ismember(extraTracks,idxMat,'rows'),:);

end

