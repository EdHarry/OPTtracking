function [coords2D, coords] = generate2DCoordsFrom3D( varargin )
%GENERATE2DCOORDSFROM3D Summary of this function goes here
%   EHarry Nov 2014

%% PARAMETERS

% defaults
param = struct('zNoise',0.01,'virtRad',0.1);

% user-selected param
var = fieldnames(param);
for i = 1:2:length(varargin)
    if ~isempty(strcmpi(varargin{i},var))
        param.(varargin{i}) = varargin{i+1};
    end
end

% unpack param
zNoise = param.zNoise;
virtRad = param.virtRad;

%% SIMULATION

% simulate 3D coords
coords = generateDataPoints_3D(varargin{:});

%% OUTPUT CONSTRUCTION

% construct 2D coords
coords2D = repmat(coords,1,2);
for i = 1:2
    j = mod(i,2) + 1;
    for iFrame = 1:size(coords2D,1)
        coords2D(iFrame,i).coords(:,j) = [];
    end
end

%% OBSCURE FEATURES

if virtRad
    for i = 1:2
        j = mod(i,2) + 1;
        for iFrame = 1:size(coords2D,1)
            disMat = full(createSparseDistanceMatrix(coords2D(iFrame,i).coords,coords2D(iFrame,i).coords,virtRad,0));
            [iRow,iCol] = find(disMat);
            idx = unique(sort([iRow iCol],2),'rows');
            for iEl = 1:size(idx,1)
                if isnan(coords2D(iFrame,j).coords(idx(iEl,1),1)) || isnan(coords2D(iFrame,j).coords(idx(iEl,2),1))
                    % do nothing
                elseif coords2D(iFrame,j).coords(idx(iEl,1),1) < coords2D(iFrame,j).coords(idx(iEl,2),1)
                    coords2D(iFrame,i).coords(idx(iEl,2),:) = NaN;
                else
                    coords2D(iFrame,i).coords(idx(iEl,1),:) = NaN;
                end
            end
        end
    end
    for i = 1:2
        for iFrame = 1:size(coords2D,1)
            coords2D(iFrame,i).coords(isnan(coords2D(iFrame,i).coords(:,1)),:) = [];
        end
    end
end

%% Z NOISE

if zNoise
    for i = 1:2
        for iFrame = 1:size(coords2D,1)
            coords2D(iFrame,i).coords(:,2) = coords2D(iFrame,i).coords(:,2) + (zNoise * randn(size(coords2D(iFrame,i).coords(:,2))));
        end
    end
end


end

