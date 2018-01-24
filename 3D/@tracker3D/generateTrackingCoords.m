function tracks = generateTrackingCoords( tracks, coords3D )
%GENERATETRACKINGCOORDS generates tracks with specified coords
%   EHarry Nov 2014

type = size(coords3D,2);
switch type
    case 1
        nDim = 3;
        tracks = placeCoords(tracks,coords3D);
    case 2
        nDim = 2;
        for i = 1:2
            tracks(i).trackInfo = placeCoords(tracks(i).trackInfo,coords3D(:,i));
        end
end

%% SUBFUNCTIONS

    function trackInfo = placeCoords(trackInfo,coords)
        for iTrack = 1:length(trackInfo)
            feats = trackInfo(iTrack).tracksFeatIndxCG;
            startTime = trackInfo(iTrack).seqOfEvents(1,1);
            idx = 0;
            for feat = feats
                idx = idx + 1;
                toModify = feat ~= 0;
                featIdx = ((idx-1)*8) + 1;
                trackInfo(iTrack).tracksCoordAmpCG(toModify,featIdx:(featIdx+nDim-1)) = coords(startTime+idx-1).coords(feat(toModify),:);
            end
        end
    end

end

