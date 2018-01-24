function output = alignFrames_3Drotation( coords3D, aor, angle, axisAngle, keepMovieInfo )
%ALIGNFRAMES aligns XYZ coordinate data based on a rotation around the axis of rotation in each frame
%aor: axis of rotation, either [1 X 2] array or [nFrames X 2] array of aor position in [x y] in each
%frame
% angle: either a single angle (radians) to rotate each frame or a
% [(nFrames-1) X 1] array of angles
% axis angle: angle between axis
%   EHarry Dec 2014

%% SETUP

if nargin < 5
    keepMovieInfo = false;
end

% get nFrames
nFrames = size(coords3D,1);
movieInfo = [];

% accept either coords3D structure or movieInfo structure
if isfield(coords3D,'xCoord')
    movieInfo = coords3D;
    movieInfo_save = movieInfo;
    isMovieInfo = true;
    coords3DmovieInfoConv('c');
else
    isMovieInfo = false;
end

%% PROJECTION ONTO ORTHONORMAL AXIS
a = pi - axisAngle;
r = [1 -cot(a); 0 csc(a)];

for iFrame = 1:nFrames
    coords3D(iFrame).coords(:,1:2) = (r * coords3D(iFrame).coords(:,1:2)')';
end

aor = (r * aor')';

%% AOR SHIFT
if size(aor,1) == 1
    aor = repmat(aor,nFrames,1);
end

for iFrame = 1:nFrames
    coords3D(iFrame).coords(:,1:2) = bsxfun(@minus,coords3D(iFrame).coords(:,1:2),aor(iFrame,:));
end

%% AOR ROTATION
if numel(angle) == 1
    angle = repmat(angle,nFrames-1,1);
end
angle = cumsum([0; angle]);
angle = reshape(angle,[1 1 nFrames]);
ca = cos(angle);
sa = sin(angle);
msa = -sa;
z = zeros(size(angle));
o = ones(size(angle));
r = [ca msa z; sa ca z; z z o];
quat = quaternion.rotationmatrixEigen(r);

for iFrame = 1:nFrames
    coords3D(iFrame).coords = quat(iFrame,:).RotateVectorCross(coords3D(iFrame).coords);
end

%% OUTPUT

if isMovieInfo && keepMovieInfo
    coords3DmovieInfoConv('m');
    for iFrame = 1:nFrames
        movieInfo(iFrame).xCoord(:,2) = movieInfo_save(iFrame).xCoord(:,2);%#ok<AGROW>
        movieInfo(iFrame).yCoord(:,2) = movieInfo_save(iFrame).yCoord(:,2);%#ok<AGROW>
        movieInfo(iFrame).zCoord(:,2) = movieInfo_save(iFrame).zCoord(:,2);%#ok<AGROW>
        movieInfo(iFrame).amp = movieInfo_save(iFrame).amp;%#ok<AGROW>
    end
    output = movieInfo;
else
    output = coords3D;
end

%% SUBFUNCTIONS

    function coords3DmovieInfoConv(desiredOutput)
        switch desiredOutput
            case 'c'
                coords3D = repmat(struct('coords',[]),nFrames,1);
                for jFrame = 1:nFrames
                    coords3D(jFrame).coords = [movieInfo(jFrame).xCoord(:,1) movieInfo(jFrame).yCoord(:,1) movieInfo(jFrame).zCoord(:,1)];
                end
                
            case 'm'
                movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),nFrames,1);
                for jFrame = 1:nFrames
                    movieInfo(jFrame).xCoord(:,1) = coords3D(jFrame).coords(:,1);
                    movieInfo(jFrame).yCoord(:,1) = coords3D(jFrame).coords(:,2);
                    movieInfo(jFrame).zCoord(:,1) = coords3D(jFrame).coords(:,3);
                end
        end
    end

end

