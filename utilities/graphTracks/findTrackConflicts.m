function conflicts = findTrackConflicts( idxMat )
%FINDTRACKCONFLICTS Summary of this function goes here
%   Detailed explanation goes here

[nTracks,nFrames] = size(idxMat);

confMat = false(nTracks);

for iFrame = 1:nFrames
    col = idxMat(:,iFrame);
    col_nZ = col(col ~= 0);
    col_unique = unique(col_nZ);
    if length(col_unique) ~= length(col_nZ)
        for idx = col_unique'
            ids = idx == col;
            if sum(ids) > 1
                ids = find(ids);
                nIds = length(ids);
                for i = 1:(nIds-1)
                    for j = (i+1):nIds
                        confMat(ids(i),ids(j)) = true;
                        confMat(ids(j),ids(i)) = true;
                    end
                end
            end
        end
    end
end

conflicts = floodFillGraph(confMat);

idx = 0;
while idx < length(conflicts)
    idx = idx + 1;
    if length(conflicts(idx).graph) == 1
        conflicts(idx) = [];
        idx = idx - 1;
    end
end

end

