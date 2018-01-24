function trackMat = constructTrackMatFromCoords( coords )
%CONSTRUCTTRACKMATFROMCOORDS Summary of this function goes here
%   Detailed explanation goes here

[nFrames,nSeries] = size(coords);

switch nSeries
    
    case 1
        try
            coords = catStruct(2,'coords.coords');
        catch
            disp('cannot track coords with missing timepoints');
            trackMat = [];
            return
        end
        
        [nTracks,nDim] = size(coords);
        nDim = nDim / nFrames;
        
        trackMat = NaN(nTracks,nFrames * 8);
        
        for iDim = 1:nDim
            trackMat(:,iDim:8:end) = coords(:,iDim:nDim:end);
        end

    case 2
        tmp = coords(:,1);%#ok<NASGU>
        tmp = catStruct(1,'tmp.coords');
        nTracks1 = size(tmp,1);
        tmp = coords(:,2);%#ok<NASGU>
        tmp = catStruct(1,'tmp.coords');
        nTracks2 = size(tmp,1);
        
        trackMat1 = sparse(nTracks1,nFrames * 8);
        trackMat2 = sparse(nTracks2,nFrames * 8);
        
        [idx1,idx2] = deal(1);
        for iFrame = 1:nFrames
            coords1 = coords(iFrame,1).coords;
            coords2 = coords(iFrame,2).coords;
            nC1 = size(coords1,1);
            nC2 = size(coords2,1);
            trackMat1(idx1:(idx1+nC1-1),((iFrame-1)*8)+1:((iFrame-1)*8)+2) = coords1;
            trackMat2(idx2:(idx2+nC2-1),((iFrame-1)*8)+1:((iFrame-1)*8)+2) = coords2;
            idx1 = idx1 + nC1;
            idx2 = idx2 + nC2;
        end
        
        trackMat = trackMat1;
        visuTracks2D(trackMat2);
        
    otherwise
        error('constructTrackMatFromCoords:unsupportedType','unsupported coord type');
end

end

