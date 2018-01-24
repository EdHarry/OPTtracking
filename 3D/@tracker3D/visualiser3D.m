function visualiser3D( tr, selectString )
%VISUALISER2D Summary of this function goes here
%   Detailed explanation goes here

filtered = ismember('f',selectString);

matlabVisVersion = ismember('m',selectString);

if ismember('2',selectString)
    
    if ismember('r',selectString)
        if filtered
            visuTracks2D(tr.tracks_filtered_raw, matlabVisVersion);
        else
            visuTracks2D(tr.tracks_raw, matlabVisVersion);
        end
    elseif ismember('t',selectString)
        if filtered
            visuTracks2D(tr.tracks_filtered_tracking, matlabVisVersion);
        else
            visuTracks2D(tr.tracks_tracking, matlabVisVersion);
        end
    else
        if filtered
            visuTracks2D(tr.tracks_filtered, matlabVisVersion);
        else
            visuTracks2D(tr.tracks, matlabVisVersion);
        end
    end
    
else
    
    if ismember('r',selectString)
        if filtered
            visuTracks3D(tr.tracks_filtered_raw, matlabVisVersion);
        else
            visuTracks3D(tr.tracks_raw, matlabVisVersion);
        end
    elseif ismember('t',selectString)
        if filtered
            visuTracks3D(tr.tracks_filtered_tracking, matlabVisVersion);
        else
            visuTracks3D(tr.tracks_tracking, matlabVisVersion);
        end
    else
        if filtered
            visuTracks3D(tr.tracks_filtered, matlabVisVersion);
        else
            visuTracks3D(tr.tracks, matlabVisVersion);
        end
    end
    
end

end

