function toRemove = filterTracks3D( coordMat, previousToRemove, direction )
%FILTERTRACKS3D Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    previousToRemove = NaN;
end

if nargin < 3
    direction = true;
end

coordMat_original = coordMat;

sigma = 6;

nTimePoints = size(coordMat,2)/8;
indexVec = false(nTimePoints,1);

vec = cat(3,coordMat(:,1:8:end),coordMat(:,2:8:end),coordMat(:,3:8:end));
vec = diff(vec,1,2);
vec = sqrt(sum(vec.^2,3));
vec2 = diff(vec,1,2);

toRemove = {};
for i = 1:size(vec2,1)
    outliers = calcOutliers(i);
    if ~isempty(outliers)
        toRemove = [toRemove; {outliers , i}];%#ok<AGROW>
        indVecTmp = indexVec;
        indVecTmp(outliers) = true;
        cc = bwconncomp(indVecTmp);
        for j = outliers'
            coordMat(i,((j-1)*8) + (1:8)) = NaN;
        end
        for j = 1:cc.NumObjects
            timepoints = cc.PixelIdxList{j};
            if min(timepoints) > 1 && max(timepoints) < nTimePoints
                n = length(timepoints);
                endTimePoint = max(timepoints) + 1;
                startTimePoint = min(timepoints) - 1;
                vector = coordMat(i,((endTimePoint-1)*8) + (1:3)) - coordMat(i,((startTimePoint-1)*8) + (1:3));
                for k = 1:n
                    coordMat(i,((startTimePoint+k-1)*8) + (1:3)) = coordMat(i,((startTimePoint-1)*8) + (1:3)) + (vector * (k / (n+1)));
                end
            end
        end
    end
end

if isequal(toRemove,previousToRemove)
    return
end

if ~isempty(toRemove)
    toRemove2 = filterTracks3D(coordMat, toRemove);
    toRemove = joinToR(toRemove,toRemove2);
end

if ~iscell(previousToRemove) && isnan(previousToRemove) && direction
    tmp = reshape(1:size(coordMat,2),8,[])';
    tmp = tmp(end:-1:1,:);
    tmp = tmp';
    tmp = tmp(:);
    toRemove2 = filterTracks3D(coordMat_original(:,tmp),NaN,false);
    if isempty(toRemove2)
        return
    elseif isempty(toRemove)
        toRemove = toRemove2;
    else
        toRemove = joinToR(toRemove, reverseToR(toRemove2), true);
    end
end



    function outliers = calcOutliers(i)
        data = vec2(i,:)';
        currentNaN = find(isnan(data));
        outliers = 0;
        while ~isempty(outliers)
            [~,~,~,outliers] = robustMean(data,[],sigma);
            data(outliers) = NaN;
        end
        outliers1 = setdiff(find(isnan(data)),currentNaN) + 1;
        data = vec(i,:)';
        currentNaN = find(isnan(data));
        outliers = 0;
        while ~isempty(outliers)
            [~,~,~,outliers] = robustMean(data,[],sigma);
            data(outliers) = NaN;
        end
        outliers = intersect(setdiff(find(isnan(data)),currentNaN),outliers1) + 1;
    end

    function toR = joinToR(toR,toR2,newTPs)
        if nargin < 3
            newTPs = false;
        end
        
        tp = cell2mat(toR(:,2));
        for ii = 1:size(toR2,1)
            jj = toR2{ii,2} == tp;
            if newTPs
                if any(jj)
                    toR{jj,1} = unique([toR{jj,1}; toR2{ii,1}]);
                else
                    toR = [toR; toR2(ii,:)];%#ok<AGROW>
                end
            elseif any(jj)
                toR{jj,1} = unique([toR{jj,1}; toR2{ii,1}]);                
            end
        end
    end

    function toR = reverseToR(toR)
        tmp = (nTimePoints:-1:1)';
        for ii = 1:size(toR,1)
            toR{ii,1} = tmp(toR{ii,1});
        end
    end

end

