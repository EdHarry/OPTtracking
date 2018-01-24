function visualiser2D( tr, selectString )
%VISUALISER2D Summary of this function goes here
%   Detailed explanation goes here

matlabVisVersion = ismember('m',selectString);

if ismember('3',selectString)
    if isempty(tr.tracks3D)
        disp('no 3D features have been generated yet');
        return
    else
        visuTracks2D(tr.tracks3D, matlabVisVersion);
    end
else
    
    if ismember('2',selectString)
        coordIdx = 2;
    else
        coordIdx = 1;
    end
    
    aligned = ismember('a',selectString);
    split = ismember('s',selectString);
    
    if aligned
        if split
            visuTracks2D(tr.split_tracks2D_aligned(coordIdx).trackInfo, matlabVisVersion);
        else
            visuTracks2D(tr.tracks2D_aligned(coordIdx).trackInfo, matlabVisVersion);
        end
    else
        if split
            visuTracks2D(tr.split_tracks2D(coordIdx).trackInfo, matlabVisVersion);
        else
            visuTracks2D(tr.tracks2D(coordIdx).trackInfo, matlabVisVersion);
        end
    end
end
end

