function coords3D = alignFrames3D_rotation( coords3D, tracks3D, verbose)
%ALIGNFRAMES aligns XY coordinate data based on estimated frame to frame displacement and rotation
% input:    cellCoords: nFrames x 2 struct with field 'coords', an nx2
%                       array, [{X,Z} , {Y,Z}]
%           tracks2D: 2D tracking information, output from 'tracking_2D'
%   EHarry Nov 2014

%% INITIALISE

if nargin < 3
    verbose = false;
end

nFrames = size(coords3D,1);

% get feature matrix
for i = 1:2
    [~,tracks] = convStruct2MatIgnoreMS(tracks3D);
    for iFrame = size(tracks,2)+1 : nFrames
        tracks = [tracks zeros(size(tracks,1),1)];%#ok<AGROW>
    end
end

%% MAIN LOOP

options = optimset('Jacobian','on','Display','off','TolX', 1e-10, 'Tolfun', 1e-10,'MaxFunEvals', 1e6, 'MaxIter', 1e6); %minimization option - use analytical Jacobian
x0 = zeros(6,1); %initial guess

for iFrame = 2:nFrames
    [coords,coordsTarget] = getCoords(iFrame,iFrame-1);
    if size(coords,1) < 6
        continue
    end
    
    % estimate shift and rotation
    x = lsqnonlin(@calcRotResiduals,x0,[],[],options,coords,coordsTarget);
    
    if verbose
        fprintf('frame: %d to %d, estimated shift = [%g %g %g], estimated rotation Euler angles = [%g %g %g]\n',iFrame,iFrame-1,x(1),x(2),x(3),x(4),x(5),x(6));
    end
    
    % add shift
    coords = bsxfun(@plus,coords3D(iFrame).coords,x(1:3)');
    
    % rotate coords
    rotQ = generateRotationQuaternions(x(4:6));
    coords = rotQ.RotateVectorCross(coords);
    coords3D(iFrame).coords = coords;
    
end

%% local functions

    function [coords,coordsTarget] = getCoords(frame,targetFrame)
        feat = tracks(:,frame);
        featTarget = tracks(:,targetFrame);
        
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
        coords = coords3D(frame).coords(feat,:);
        coordsTarget = coords3D(targetFrame).coords(featTarget,:);
        
        % save coords
        coords_ = coords;
        coordsTarget_ = coordsTarget;
        
        % remove outliers
        % first shift by com
        coords = bsxfun(@minus,coords,mean(coords));
        coordsTarget = bsxfun(@minus,coordsTarget,mean(coordsTarget));
        meanX = robustMean(coords(:,1));
        meanY = robustMean(coords(:,2));
        meanZ = robustMean(coords(:,3));
        coords = bsxfun(@minus,coords,[meanX meanY meanZ]);
        meanX = robustMean(coordsTarget(:,1));
        meanY = robustMean(coordsTarget(:,2));
        meanZ = robustMean(coordsTarget(:,3));
        coordsTarget = bsxfun(@minus,coordsTarget,[meanX meanY meanZ]);
        
        % detect outliers
        outliers = 0;
        while ~isempty(outliers)
            sep = coordsTarget - coords;
            sep = sep';
            [~,~,~,outliers] = robustMean(sep(:));
            outliers = ceil(outliers/3);
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


    function [rotQuat,devRotSets] = generateRotationQuaternions(angles) % generates rotation quaternion (in 2D) and its derivative (relative to the angle) (note derivative(3,3) = 1 here instead of 0 since it will be used on 2D vectors)
        c = cos(angles);
        s = sin(angles);
        ms = -s;
        R = zeros([3 3 3]);
        R(:,:,1) = [c(1) ms(1) 0; s(1) c(1) 0; 0 0 1];
        R(:,:,2) = [c(2) 0 s(2); 0 1 0; ms(2) 0 c(2)];
        R(:,:,3) = [1 0 0; 0 c(3) ms(3); 0 s(3) c(3)];
        rotQuat = quaternion.rotationmatrixEigen(R(:,:,1)*R(:,:,2)*R(:,:,3));
        if nargout > 1
            mc = -c;
            Rd = zeros([3 3 3]);
            Rd(:,:,1) = [ms(1) mc(1) 0; c(1) ms(1) 0; 0 0 0];
            Rd(:,:,2) = [ms(2) 0 c(2); 0 0 0; mc(2) 0 ms(2)];
            Rd(:,:,3) = [0 0 0; 0 ms(3) mc(3); 0 c(3) ms(3)];
            devRotSets = repmat(struct('quat',[],'mat',[]),3,1);
            for index = 1:3
                otherIdx = setdiff(1:3,index);
                devRotSets(index).quat = quaternion.rotationmatrixEigen(R(:,:,otherIdx));
                devRotSets(index).mat = Rd(:,:,index);
            end
        end
    end

    function [F,J] = calcRotResiduals(x,coords,coordsTarget) % objective function
        % zero array for Jacobian
        z = zeros(size(coords,1),3);
        
        % generate rotations
        [rotQuat,devRotSets] = generateRotationQuaternions(x(4:6));
        
        % add shift to coords
        coords = bsxfun(@plus,coords,x(1:3)');
        
        % get residuals
        F = rotQuat.RotateVectorCross(coords) - coordsTarget;
        F = F(:);
        
        %calculate the Jacobian matrix - initialise
        J = zeros(length(F),6);
        
        % derivative relative to x shift
        v = z;
        v(:,1) = 1;
        jCol = rotQuat.RotateVectorCross(v);
        J(:,1) = jCol(:);
        
        % derivative relative to y shift
        v = z;
        v(:,2) = 1;
        jCol = rotQuat.RotateVectorCross(v);
        J(:,2) = jCol(:);
        
        % derivative relative to z shift
        v = z;
        v(:,3) = 1;
        jCol = rotQuat.RotateVectorCross(v);
        J(:,3) = jCol(:);
        
        % derivative relative to x Euler angle
        jCol = (devRotSets(1).mat * devRotSets(1).quat(1,:).RotateVectorCross(devRotSets(1).quat(2,:).RotateVectorCross(coords))')';
        J(:,4) = jCol(:);
        
        % derivative relative to y Euler angle
        jCol = devRotSets(2).quat(1,:).RotateVectorCross((devRotSets(2).mat * devRotSets(2).quat(2,:).RotateVectorCross(coords)')');
        J(:,5) = jCol(:);
        
        % derivative relative to z Euler angle
        jCol = devRotSets(3).quat(1,:).RotateVectorCross(devRotSets(3).quat(2,:).RotateVectorCross((devRotSets(3).mat * coords')'));  
        J(:,6) = jCol(:);
        
    end

end

