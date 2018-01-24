function trackGraphInfo = compound2Graph(seqOfEvents)
%COMPOUND2GRAPH Summary of this function goes here
%   EHarry Nov 2014

n = size(seqOfEvents,1);
trackGraph = false(n);
subTrackIdx = cell(n,1);
[startNodes,endNodes,splitNodes,mergeNodes] = deal([]);

for iNode = 1:n
    myTime = seqOfEvents(iNode,1);
    typeOfEvent = seqOfEvents(iNode,2);
    myIdx = seqOfEvents(iNode,3);
    targetIdx = seqOfEvents(iNode,4);
    
    subTrackIdx{iNode} = [subTrackIdx{iNode} myIdx];
    
    if isnan(targetIdx)
        if typeOfEvent == 1
            startNodes = [startNodes; iNode];%#ok<AGROW>
        else
            endNodes = [endNodes; iNode];%#ok<AGROW>
        end
    else
        if typeOfEvent == 1
            splitNodes = [splitNodes; iNode];%#ok<AGROW>
            idx = 0;
            loop = true;
            while loop
                idx = idx - 1;
                jNode = iNode + idx;
                if ismember(targetIdx,subTrackIdx{jNode}) && myTime ~= seqOfEvents(jNode,1)
                    loop = false;
                    subTrackIdx{jNode} = [subTrackIdx{jNode} myIdx];
                    trackGraph(jNode,iNode) = true;
                end
            end
        else
            mergeNodes = [mergeNodes; iNode];%#ok<AGROW>
            subTrackIdx{iNode} = [subTrackIdx{iNode} targetIdx];
            idx = 0;
            loop = true;
            while loop
                idx = idx - 1;
                jNode = iNode + idx;
                if ismember(myIdx,subTrackIdx{jNode}) && myTime ~= seqOfEvents(jNode,1)
                    loop = false;
                    trackGraph(jNode,iNode) = true;
                end
            end
        end
    end
end

uniqueSubIdxs = unique(seqOfEvents(:,3))';
for subIdx = uniqueSubIdxs
    idx = 0;
    loop = true;
    while loop
        idx = idx + 1;
        if ismember(subIdx,subTrackIdx{idx})
            loop = false;
        end
    end
    
    loop = true;
    while loop
        jdx = idx;
        loop2 = true;
        while loop2
            jdx = jdx + 1;
            if jdx > n || ismember(subIdx,subTrackIdx{jdx})
                loop2 = false;
            end
        end
        if jdx > n
            loop = false;
        else
            trackGraph(idx,jdx) = true;
            idx = jdx;
        end
    end
end

trackGraphInfo.trackGraph = trackGraph;
trackGraphInfo.subTrackIdx = subTrackIdx;
trackGraphInfo.startNodes = startNodes;
trackGraphInfo.endNodes = endNodes;
trackGraphInfo.splitNodes = splitNodes;
trackGraphInfo.mergeNodes = mergeNodes;

end

