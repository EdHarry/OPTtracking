function cellCoords = alignFrames_rotation( cellCoords, tracks2D, verbose, toAlign)
%ALIGNFRAMES aligns XY coordinate data based on estimated frame to frame displacement and rotation
% input:    cellCoords: nFrames x 2 struct with field 'coords', an nx2
%                       array, [{X,Z} , {Y,Z}]
%           tracks2D: 2D tracking information, output from 'tracking_2D'
%   EHarry Nov 2014

%% INITIALISE

if nargin < 3
    verbose = false;
end

if nargin < 4
    toAlign = 1:2;
end

nFrames = size(cellCoords,1);

% get feature matrix
tracksMat = repmat(struct('tracks',[]),1,2);
for i = toAlign
    [~,tracksMat(i).tracks] = convStruct2MatIgnoreMS(tracks2D(i).trackInfo);
    for iFrame = size(tracksMat(i).tracks,2)+1 : nFrames
        tracksMat(i).tracks = [tracksMat(i).tracks zeros(size(tracksMat(i).tracks,1),1)];
    end
end

%% MAIN LOOP

options = optimset('Jacobian','on','Display','off','TolX', 1e-10, 'Tolfun', 1e-10,'MaxFunEvals', 1e6, 'MaxIter', 1e6); %minimization option - use analytical Jacobian
x0 = [0; 0; 0]; %initial guess

for i = toAlign
    for iFrame = 2:nFrames
        [coords,coordsTarget] = getCoords(cellCoords(:,i),tracksMat(i).tracks,iFrame,iFrame-1);
        if size(coords,1) < 3
            continue
        end
        
        % estimate shift and rotation
        x = lsqnonlin(@calcRotResiduals,x0,[],[],options,coords,coordsTarget);
        
        if verbose
            fprintf('series: %d, frame: %d to %d, estimated shift = [%g %g], estimated rotation angle = %g\n',i,iFrame,iFrame-1,x(1),x(2),x(3));
        end
        
        % add shift
        coords = bsxfun(@plus,cellCoords(iFrame,i).coords,x(1:2)');
        
        % rotate coords
        rotQ = generateRotationQuaternions(x(3));
        coords = rotQ.RotateVectorCross([coords zeros(size(coords,1),1)]);
        cellCoords(iFrame,i).coords = coords(:,1:2);
        
    end
end

%% local functions

    function [coords,coordsTarget] = getCoords(cellCoords,trackIdx,frame,targetFrame)
        feat = trackIdx(:,frame);
        featTarget = trackIdx(:,targetFrame);
        
        %find those features that are linked between the two frames
        goodIndx = feat ~= 0 & featTarget ~= 0;
        feat = feat(goodIndx);
        featTarget = featTarget(goodIndx);
        
        % return if either featIdx is completely empty
        if all(feat==0) || all(featTarget==0)
            coords = [];
            coordsTarget = [];
            return
        end
        
        % get coords
        coords = cellCoords(frame).coords(feat,:);
        coordsTarget = cellCoords(targetFrame).coords(featTarget,:);
        
        % save coords
        coords_ = coords;
        coordsTarget_ = coordsTarget;
        
        % remove outliers
        % first shift by com
        coords = bsxfun(@minus,coords,mean(coords));
        coordsTarget = bsxfun(@minus,coordsTarget,mean(coordsTarget));
        meanX = robustMean(coords(:,1));
        meanY = robustMean(coords(:,2));
        coords = bsxfun(@minus,coords,[meanX meanY]);
        meanX = robustMean(coordsTarget(:,1));
        meanY = robustMean(coordsTarget(:,2));
        coordsTarget = bsxfun(@minus,coordsTarget,[meanX meanY]);
        
        % detect outliers
        outliers = 0;
        while ~isempty(outliers)
            sep = coordsTarget - coords;
            sep = sep';
            [~,~,~,outliers] = robustMean(sep(:));
            outliers = ceil(outliers/2);
            outliers = unique(outliers);
            coords(outliers,:) = NaN;
            coordsTarget(outliers,:) = NaN;
        end
        
        % get original coords and remove outliers
        outliers = isnan(coords(:,1));
        coords = coords_;
        coordsTarget = coordsTarget_;
        coords(outliers,:) = [];
        coordsTarget(outliers,:) = [];
    end


    function [rotQuat,devRotQuat] = generateRotationQuaternions(angle) % generates rotation quaternion (in 2D) and its derivative (relative to the angle) (note derivative(3,3) = 1 here instead of 0 since it will be used on 2D vectors)
        c = cos(angle);
        s = sin(angle);
        ms = -s;
        R = [c ms 0; s c 0; 0 0 1];
        rotQuat = quaternion.rotationmatrixEigen(R);
        if nargout > 1
            mc = -c;
            R = [ms mc 0; c ms 0; 0 0 1];
            devRotQuat = quaternion.rotationmatrixEigen(R);
        end
    end

    function [F,J] = calcRotResiduals(x,coords,coordsTarget) % objective function
        % zero array, for jacobian
        z = zeros(size(coords,1),3);
        
        % generate rotations
        [rotQuat,devRotQuat] = generateRotationQuaternions(x(3));
        
        % add shift to coords
        coords = bsxfun(@plus,coords,x(1:2)');
        
        % get residuals
        F = rotQuat.RotateVectorCross([coords z(:,1)]) - [coordsTarget z(:,1)];
        F = F(:,1:2);
        F = F(:);
        
        %calculate the Jacobian matrix - initialise
        J = zeros(length(F),3);
        
        % derivative relative to x shift
        v = z;
        v(:,1) = 1;
        jCol = rotQuat.RotateVectorCross(v);
        jCol = jCol(:,1:2);
        J(:,1) = jCol(:);
        
        % derivative relative to y shift
        v = z;
        v(:,2) = 1;
        jCol = rotQuat.RotateVectorCross(v);
        jCol = jCol(:,1:2);
        J(:,2) = jCol(:);
        
        % derivative relative to angle
        jCol = devRotQuat.RotateVectorCross([coords z(:,1)]);
        jCol = jCol(:,1:2);
        J(:,3) = jCol(:);
        
    end

end

