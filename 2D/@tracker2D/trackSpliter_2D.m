function tracks2D = trackSpliter_2D( tracks2D )
%TRACKSPLITER2D Summary of this function goes here
%   Detailed explanation goes here

for i = 1:2
    tracks2D(i).trackInfo = splitTracks2D(tracks2D(i).trackInfo);
end


end

