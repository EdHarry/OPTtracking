function [ coords3D, indexes2D ] = get3DProjections( coords2D, notToProcess )
%GET3DPROJECTIONS Summary of this function goes here
%   Detailed explanation goes here

%% SETUP
if nargin < 2
    notToProcess = [0 0];
end

%% GET PAIRING PARAMETERS
pairingParam = trackPairing_parameters;

% unpack
maxZDifference = pairingParam.maxZDifference;

%% PROCESSING
nCoords1 = size(coords2D(1).coords,1);
nCoords2 = size(coords2D(2).coords,1);

indexes2D = [];
coords3D = [];

for iCoord1 = 1:nCoords1 
    for iCoord2 = 1:nCoords2
        if ~ismember([iCoord1 iCoord2],notToProcess,'rows') && abs(coords2D(1).coords(iCoord1,2) - coords2D(2).coords(iCoord2,2)) < maxZDifference
            indexes2D = [indexes2D; iCoord1 iCoord2];%#ok<AGROW>
            coords3D = [coords3D; coords2D(1).coords(iCoord1,1) coords2D(2).coords(iCoord2,1) mean([coords2D(1).coords(iCoord1,2) coords2D(2).coords(iCoord2,2)])];%#ok<AGROW>
        end
    end
end

end

