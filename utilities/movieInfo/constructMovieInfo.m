function [ movieInfo, indexMap ] = constructMovieInfo( indexMat, coordMat, probDim )
%CONSTRUCTMOVIEINFO Summary of this function goes here
%   Detailed explanation goes here

nFrames = size(indexMat,2);
switch probDim
    case 2
        movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),nFrames,1);
    case 3
        movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),nFrames,1);
end
indexMap = repmat(struct('map',[]),nFrames,1);

for iFrame = 1:nFrames
    indexes = indexMat(:,iFrame);
    coords = coordMat(:,(8*(iFrame-1))+1:(8*(iFrame-1))+8);
    x = [coords(:,1:8:end) coords(:,5:8:end)];
    y = [coords(:,2:8:end) coords(:,6:8:end)];
    z = [coords(:,3:8:end) coords(:,7:8:end)];
    amp = [coords(:,4:8:end) coords(:,8:8:end)];
    goodCoord = ~isnan(x(:,1));
    x(isnan(x(:,2)),2) = 0;
    y(isnan(y(:,2)),2) = 0;
    z(isnan(z(:,2)),2) = 0;
    amp(isnan(amp(:,1)),:) = 0;
    movieInfo(iFrame).xCoord = x(goodCoord,:);
    movieInfo(iFrame).yCoord = y(goodCoord,:);
    if probDim == 3
        movieInfo(iFrame).zCoord = z(goodCoord,:);
    end
    movieInfo(iFrame).amp = amp(goodCoord,:);
    indexMap(iFrame).map = indexes(goodCoord);
end

end

