function toRemove = filterTracksDis3D( coordMat )
%FILTERTRACKSDIS3D Summary of this function goes here
%   Detailed explanation goes here

sigma = 6;

trackSEL = getTrackSEL(coordMat);
[trackSEL, trackSEL2] = deal(trackSEL(:,3));

outliers = 0;
while ~isempty(outliers)
    [rMean,~,~,outliers] = robustMean(trackSEL,[],sigma);
    trackSEL(outliers) = NaN;
end

toRemove = find(isnan(trackSEL) & trackSEL2 < rMean);

end

