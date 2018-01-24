function [ coords3D, indexes2D ] = findBest2DProjections( coords2D, notToProcess )
%GET3DPROJECTIONS Summary of this function goes here
%   Detailed explanation goes here

%% SETUP
if nargin < 2
    notToProcess = cell(1,2);
end

%% GET PAIRING PARAMETERS
pairingParam = trackPairing_parameters;

% unpack
maxZDifference = pairingParam.maxZDifference;

%% PROCESSING
nCoords1 = size(coords2D(1).coords,1);
nCoords2 = size(coords2D(2).coords,1);

toProcess = cell(2);
toProcess{1,2} = 1:nCoords1;
toProcess{1,1} = setdiff(1:nCoords1,notToProcess{1});
toProcess{2,2} = 1:nCoords2;
toProcess{2,1} = setdiff(1:nCoords2,notToProcess{2});

indexes2D = [];
coords3D = [];

for i = 1:2
    j = mod(i,2) + 1;
    for iCoord1 = toProcess{i,1}
        indexes2D_ind = [];
        best = Inf;
        for iCoord2 = toProcess{j,2}
            indexes2D_ind_tmp = [iCoord1 iCoord2];
            indexes2D_ind_tmp = indexes2D_ind_tmp([i j]);
            idx1 = indexes2D_ind_tmp(1);
            idx2 = indexes2D_ind_tmp(2);
            if all(abs(coords2D(1).coords(idx1,2) - coords2D(2).coords(idx2,2)) < [maxZDifference best])
                best = abs(coords2D(1).coords(idx1,2) - coords2D(2).coords(idx2,2));
                indexes2D_ind = indexes2D_ind_tmp;
            end
        end
        if ~isempty(indexes2D_ind)
            idx1 = indexes2D_ind(1);
            idx2 = indexes2D_ind(2);
            indexes2D = [indexes2D; indexes2D_ind];%#ok<AGROW>
            coords3D = [coords3D; coords2D(1).coords(idx1,1) coords2D(2).coords(idx2,1) mean([coords2D(1).coords(idx1,2) coords2D(2).coords(idx2,2)])];%#ok<AGROW>
            toProcess{j,1} = setdiff(toProcess{j,1},indexes2D_ind(j));
        end
    end
end

end

