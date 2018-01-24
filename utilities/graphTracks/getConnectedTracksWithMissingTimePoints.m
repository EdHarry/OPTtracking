function islands = getConnectedTracksWithMissingTimePoints(idxMat)
%GETCONNECTEDTRACKSWITHMISSINGTIMEPOINTS Summary of this function goes here
%   Detailed explanation goes here

[tracksToProcess,framesToProcess] = find(idxMat == 0);
nNodes = length(tracksToProcess);

if nNodes == 0
    islands = [];
    return
end

graph = false(nNodes);
for iNode = 1:(nNodes-1)
    myTrackId = tracksToProcess(iNode);
    myFrameId = framesToProcess(iNode);
    jNode = iNode;
    loop = true;
    while jNode < nNodes && loop
        jNode = jNode + 1;
        if myTrackId == tracksToProcess(jNode)
            [graph(iNode,jNode),graph(jNode,iNode)] = deal(true);
            loop = false;
        end
    end
    jNode = iNode;
    loop = true;
    while jNode < nNodes && loop
        jNode = jNode + 1;
        if myFrameId == framesToProcess(jNode)
            [graph(iNode,jNode),graph(jNode,iNode)] = deal(true);
            loop = false;
        end
    end
end
islands = floodFillGraph(graph);

for iIsland = 1:length(islands)
    islands(iIsland).graph = unique(tracksToProcess(islands(iIsland).graph))';
end

end

