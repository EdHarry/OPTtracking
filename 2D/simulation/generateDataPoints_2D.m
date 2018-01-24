function coords = generateDataPoints_2D( varargin )
%GENERATEDATAPOINTS_2D Summary of this function goes here
%   Detailed explanation goes here

%% SETUP

% default param
param = struct('nFrames',100,'nFeat',10,'initCoord',5,'initVel',0.1,'velVar',0.05,'angVar',0.05,'radius',5,'repulseF',0.1,'mergeProb',0.001,'splitProb',0.001);

% user-selected param
var = fieldnames(param);
for i = 1:2:length(varargin)
    if ~isempty(strcmpi(varargin{i},var))
        param.(varargin{i}) = varargin{i+1};
    end
end

% unpack param
nFrames = param.nFrames;
nFeat = param.nFeat;
initCoord = param.initCoord;
initVel = param.initVel;
velVar = param.velVar;
angVar = param.angVar;
radius = param.radius;
repulseF = param.repulseF;
mergeProb = param.mergeProb;
splitProb = param.splitProb;

%% INITIATE

if nFrames < 3
    error('generateDataPoints:nFrames','error, nFrames must be > 2');
end

if isempty(initCoord) % generate random starting locations
    initCoord = randn(nFeat,2);
    initCoord = bsxfun(@minus,initCoord,mean(initCoord));
elseif numel(initCoord) == 1 && isnumeric(initCoord)
    initCoord = initCoord * randn(nFeat,2);
    initCoord = bsxfun(@minus,initCoord,mean(initCoord));
else
    if ~isnumeric(initCoord) || size(initCoord,2) ~= 2
        error('generateDataPoints:initCoord','initCoord not of correct type');
    end
    nFeat = size(initCoord,1);
end

if isempty(initVel)
    initVel = randn(nFeat,2);
    [~,initVel] = normList(initVel);
elseif numel(initVel) == 1 && isnumeric(initVel)
    tmp = randn(nFeat,2);
    [~,tmp] = normList(tmp);
    initVel = initVel * tmp;
    clear tmp
elseif any(size(initVel) ~= [nFeat 2]) || ~isnumeric(initVel)
    error('generateDataPoints:initVel','initVel not of correct type');
end

coords = repmat(struct('coords',[]),nFrames,1);
coords(1).coords = initCoord;

currentV = initVel;

[repulsions,acc] = deal(zeros(nFeat,2));

splitIdx = repmat(struct('idx',[]),nFrames,1);
tmpPrevCoords = [];

%% SIMULATE

% velocity of first point
iFrame = 2;
coords(iFrame).coords = coords(iFrame-1).coords + currentV;
previousV = currentV;
updateVelocity;
calcRepultions;

% loop
for iFrame = 3:nFrames
    acc = currentV - previousV + repulsions; % calc accelerations
    coords(iFrame).coords = acc + (2 * coords(iFrame-1).coords) - [coords(iFrame-2).coords; tmpPrevCoords];
    split;
    updateVelocity;
    calcRepultions;
end

% merging
if mergeProb
    for iFrame = nFrames:-1:3
        merge;
    end
end

%% SUBFUNCIONS

    function updateVelocity
        [mag,tmp] = normList(previousV + acc);
        previousV = currentV;
        currentV = tmp;
        mag = mag + (velVar * randn(size(currentV,1),1));
        currentV = bsxfun(@times,currentV,mag);
        currentV = generateRandomRotation(currentV);
    end

    function calcRepultions
        repulsions = zeros(nFeat,2);
        disMat = full(createSparseDistanceMatrix(coords(iFrame).coords,coords(iFrame).coords,radius,0));
        [iRow,iCol] = find(disMat);
        idxs = unique(sort([iRow iCol],2),'rows');
        iRow = idxs(:,1);
        iCol = idxs(:,2);
        for iEl = 1:size(idxs,1)
            [~,dirVec] = normList(coords(iFrame).coords(iRow(iEl),:) - coords(iFrame).coords(iCol(iEl),:));
            repMag = ((radius - disMat(iRow(iEl),iCol(iEl))) / radius) * repulseF;
            repulsions(iRow(iEl),:) = repulsions(iRow(iEl),:) + (repMag * dirVec);
            repulsions(iCol(iEl),:) = repulsions(iCol(iEl),:) - (repMag * dirVec);
        end
    end

    function split
        if splitProb
            toSplit = splitProb > rand(nFeat,1);
            if any(toSplit)
                coordsToSplit_m1 = coords(iFrame-1).coords(toSplit,:);
                tmp = [coords(iFrame-2).coords; tmpPrevCoords];
                coordsToSplit_m2 = tmp(toSplit,:);
                repToSplit = repulsions(toSplit,:);
                preVelToSplit = previousV(toSplit,:);
                curVelToSplit = currentV(toSplit,:);
                %accToSplit = acc(toSplit,:);
                curVelToSplit = generateRandomRotation(curVelToSplit);
                accToSplit = curVelToSplit - preVelToSplit + repToSplit;
                %accToSplit = generateRandomRotation(accToSplit);
                splitCoords = accToSplit + (2 * coordsToSplit_m1) - coordsToSplit_m2;
                %curVelToSplit = accToSplit + preVelToSplit - repToSplit;
                splitIdx(iFrame).idx = (size(coords(iFrame).coords,1) + 1) : (size(splitCoords,1) + size(coords(iFrame).coords,1));
                coords(iFrame).coords = [coords(iFrame).coords; splitCoords];
                currentV = [currentV; curVelToSplit];
                previousV = [previousV; preVelToSplit];
                acc = [acc; accToSplit];
                nFeat = size(currentV,1);
                tmpPrevCoords = coordsToSplit_m1;
            else
                tmpPrevCoords = [];
            end
        end
    end

    function merge
        toMerge = mergeProb > rand(size(coords(iFrame).coords,1),1);
        toMerge(splitIdx(iFrame).idx) = false;
        toMergeIdx = find(toMerge);
        if any(toMerge)
            newCoord = coords(iFrame).coords(toMerge,:);
            newCoord_all = repmat(struct('coords',[]),iFrame-1,1);
            idx = 0;
            for jFrame = iFrame:-1:3
                idx = idx + 1;
                toMerge(splitIdx(jFrame).idx) = false;
                coordsToMerge = newCoord(ismember(toMergeIdx,find(toMerge)),:);
                coordsToMerge_m1 = coords(jFrame-1).coords(toMerge,:);
                coordsToMerge_m2 = coords(jFrame-2).coords(toMerge,:);
                accToMerge = coordsToMerge - (2 * coordsToMerge_m1) + coordsToMerge_m2;
                if jFrame == iFrame
                    accToMerge = generateRandomRotation(accToMerge);
                end
                newCoord = (coordsToMerge + coordsToMerge_m2 - accToMerge) / 2;
                newCoord_all(idx).coords = newCoord;
            end
            vel = coords(2).coords(toMerge,:) - coords(1).coords(toMerge,:);
            newCoord_all(idx+1).coords = newCoord - vel;
            for jFrame = 1:(iFrame-1)
                kFrame = iFrame - jFrame;
                coords(kFrame).coords = [coords(kFrame).coords; newCoord_all(jFrame).coords];
            end
        end
    end

    function vec = generateRandomRotation(vec)
        n = size(vec,1);
        angles = angVar * randn(n,1) * 2 * pi;
        angles = reshape(angles,1,1,size(angles,1));
        c = cos(angles);
        s = sin(angles);
        ms = -s;
        z = zeros(size(angles));
        o = ones(size(angles));
        r = [c ms z; s c z; z z o];
        quat = quaternion.rotationmatrixEigen(r);
        vec = quat.RotateVectorCross([vec zeros(n,1)]);
        vec = vec(:,1:2);
    end

end

