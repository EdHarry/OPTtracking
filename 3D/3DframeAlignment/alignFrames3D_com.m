function coords3D = alignFrames3D_com( coords3D )
%ALIGNFRAMES_COM align cellCoords based on a linear shift around the centre
%of mass in each frame
%   EHarry Nov 2014

nFrames = size(coords3D,1);

for iFrame = 1:nFrames
    coords = coords3D(iFrame).coords;
    com = mean(coords);
    coords3D(iFrame).coords = bsxfun(@minus,coords,com);
end

end

