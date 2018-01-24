function [pathsOutput, earlyExit] = findAllPaths( graphMat, startIdx, endIdx )
%FINDALLPATHS Summary of this function goes here
%   Detailed explanation goes here

%% SETUP

maxGraphSize = 10;
maximumMatrixProcessingStats = [0.2 56];
earlyExit = false;

nNodes = size(graphMat,1);
if mod(nNodes,2) == 1
    graphMat = [graphMat false(nNodes,1) ; false(1,nNodes+1)];
    nNodes = nNodes + 1;
end
nNodes = nNodes / 2;
[subStartIdxs,subEndIdxs] = deal((1:nNodes)');

%% EXIT IF MATRIX VIOLATES MAXIMUM PROCESSING STATS

if nNodes >= maximumMatrixProcessingStats(2)/2 && nnz(graphMat)/numel(graphMat) >= maximumMatrixProcessingStats(1)
    pathsOutput = {};
    earlyExit = true;
    return
end

%% PROCESS SUBGRAPHS

paths = cell(2);
for i = 1:2
    for j = 1:2
        iIndex = (((i-1) * nNodes) + 1) : (((i-1) * nNodes) + nNodes);
        jIndex = (((j-1) * nNodes) + 1) : (((j-1) * nNodes) + nNodes);
        subGraphMat = graphMat(iIndex,jIndex);
        if i == j
            if nNodes > maxGraphSize
                [extraPaths, earlyExit] = findAllPaths(subGraphMat, subStartIdxs, subEndIdxs);
                if earlyExit
                    pathsOutput = {};
                    earlyExit = true;
                    return
                end
            else
                N = nNodes;
                extraPaths = getPaths;
            end
            if i == 2
                extraPaths = cellfun(@(x)(x + nNodes),extraPaths,'UniformOutput',false);
            end
        else
            N = 2;
            extraPaths = getPaths;
            diagPaths = find(diag(subGraphMat));
            if ~isempty(diagPaths)
                extraPaths = [extraPaths; mat2cell([diagPaths diagPaths],ones(length(diagPaths),1),2)];%#ok<AGROW>
            end
            switch j
                case 1
                    extraPaths = cellfun(@(x)([x(:,1)+nNodes x(:,2)]),extraPaths,'UniformOutput',false);
                case 2
                    extraPaths = cellfun(@(x)([x(:,1) x(:,2)+nNodes]),extraPaths,'UniformOutput',false);
            end
        end
        paths{i,j} = extraPaths;
    end
end

%% EXIT IF NO PATHS FOUND

if all(cellfun(@isempty,paths(:)))
    pathsOutput = {};
    return
end

%% SUBMAT NODE TRANSITION PATHS

subMatTransitionPaths = cell(4);
paths = paths';
N = 2;
for i = 1:4
    for j = 1:4
        if i ~= j
            endIdx1 = cellfun(@(x)(x(end)),paths{i});
            startIdx2 = cellfun(@(x)(x(1)),paths{j});
            nEnd = length(endIdx1);
            nStart = length(startIdx2);
            subNnodes = nEnd + nStart;
            [idxAll,idAll] = deal(zeros(min(100000000,nEnd*nStart),1));
            k = 0;
            for idx = 1:nEnd
                k = k + 1;
                id = find(ismember(startIdx2,endIdx1(idx))) + nEnd;
                nId = length(id);
                idxId = k:(nId+k-1);
                idxAll(idxId) = repmat(idx,nId,1);
                idAll(idxId) = id;
                k = nId+k-1;
            end
            idxId = idxAll ~= 0;
            idxAll = idxAll(idxId);
            idAll = idAll(idxId);
            clear idxId
            subGraphMat = sparse(idxAll,idAll,true(length(idxAll),1),subNnodes,subNnodes);
            subStartIdxs = (1:nEnd)';
            subEndIdxs = ((nEnd+1):subNnodes)';
            subPaths = cell2mat(getPaths);
            clear subGraphMat idAll idxAll
            if ~isempty(subPaths)
                subPaths(:,2) = subPaths(:,2) - nEnd;
            else
                subPaths = [];
            end
            subMatTransitionPaths{i,j} = subPaths;
        end
    end
end

%% CYCLIC SUBMAT TRANSITIONS

subMat = ~cellfun(@isempty,subMatTransitionPaths);

% single (no) transitions
pathsOutput = {};
currentTrackList = {};
currentSelect = {};
subTrackIds = [];
for i = 1:4
    subPaths = paths{i};
    select = find(filterByStart);
    subPaths = subPaths(select);
    if ~isempty(subPaths)
        pathsOutput = [pathsOutput; subPaths(filterByEnd)];%#ok<AGROW>
        currentTrackList = [currentTrackList; {subPaths}];%#ok<AGROW>
        subTrackIds = [subTrackIds; i];%#ok<AGROW>
        currentSelect = [currentSelect; {select}];%#ok<AGROW>
    end
end

% exit if no possible starts
if isempty(subTrackIds)
    pathsOutput = {};
    return
end

N = 1;
subEndIdxs = (1:4)';
subStartIdxs = subTrackIds;
loop = true;
while loop
    N = N + 1;
    subStartIdxs = cell2mat(allPaths(subMat, subStartIdxs, subEndIdxs, N, N, true));
    pathAddedThisN = false;
    idx = 0;
    [newSelect,newPaths] = deal(cell(size(subStartIdxs,1),1));
    while idx < size(subStartIdxs,1)
        idx = idx + 1;
        path = subStartIdxs(idx,:);
        node = path(end);
        previousNode = path(end-1);
        selectCurrentTracks = ismember(subTrackIds,path(1:end-1),'rows');
        subPaths = currentTrackList{selectCurrentTracks};
        select = currentSelect{selectCurrentTracks};
        transition_ = subMatTransitionPaths{previousNode,node};
        nSelect = length(select);
        transition = zeros(min(100000000,nSelect*size(transition_,1)),2);
        id = 0;
        for k = 1:nSelect
            id = id + 1;
            tmp = transition_(transition_(:,1)==select(k),2);
            nTmp = size(tmp,1);
            tmp = [repmat(k,nTmp,1) tmp];%#ok<AGROW>
            transition(id:(nTmp+id-1),:) = tmp;
            id = nTmp+id-1;
        end
        transition = transition(transition(:,1)~=0,:);
        select = transition(:,2);
        tmp = paths{node}(transition(:,2));
        subPaths = cellfun(@(x)(x(1:end-1)),subPaths(transition(:,1)),'UniformOutput',false);
        subPaths = cellfun(@(x,y)([x y]),subPaths,tmp,'UniformOutput',false);
        goodPaths = cellfun(@(x)(length(unique(x))==length(x)),subPaths);
        subPaths = subPaths(goodPaths);
        select = select(goodPaths);
        newSelect{idx} = select;
        newPaths{idx} = subPaths;
        if isempty(subPaths)
            subStartIdxs(idx,:) = [];
            newSelect(idx) = [];
            newPaths(idx) = [];
            idx = idx - 1;
        else
            pathAddedThisN = true;
            pathsOutput = [pathsOutput; subPaths(filterByEnd)];%#ok<AGROW>
        end
    end
    if ~pathAddedThisN
        loop = false;
    else
        subTrackIds = subStartIdxs;
        currentTrackList = newPaths;
        currentSelect = newSelect;
    end
end

if isempty(pathsOutput)
    pathsOutput = {};
end

%% SUBFUNCTIONS

    function paths = getPaths
        paths = {};
        p = allPaths(subGraphMat, subStartIdxs, subEndIdxs, N);
        for i_ = 1:size(p,1)
            p_ = p{i_,1};
            paths = [paths; mat2cell(p_,ones(size(p_,1),1),size(p_,2))];%#ok<AGROW>
        end
    end

    function select = filterByStart
        select = cellfun(@(x)(ismember(x(1),startIdx)),subPaths);
    end

    function select = filterByEnd
        select = cellfun(@(x)(ismember(x(end),endIdx)),subPaths);
    end

end

