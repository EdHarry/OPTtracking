function [cellCoords,coms] = alignFrames_com( cellCoords )
%ALIGNFRAMES_COM align cellCoords based on a linear shift around the centre
%of mass in each frame
%   EHarry Nov 2014

nFrames = size(cellCoords,1);
coms = repmat(struct('com',zeros(1,2)),nFrames,2);

for iFrame = 1:nFrames
    for i = 1:2
        coords = cellCoords(iFrame,i).coords;
        com = mean(coords);
        cellCoords(iFrame,i).coords = bsxfun(@minus,coords,com);
        coms(iFrame,i).com = com;
    end
end

end

