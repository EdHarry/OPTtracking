
%% general gap closing parameters
%
% Copyright (C) 2014 LCCB
%
% This file is part of u-track.
%
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
%
gapCloseParam.timeWindow = 7; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%parameters

parameters.linearMotion = 2; %use linear motion Kalman filter.

parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 4.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = 2; %use linear motion Kalman filter.

parameters.minSearchRadius = 2; %minimum allowed search radius.
parameters.maxSearchRadius = 4.5; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0.25 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

%parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
%saveResults.dir = 'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\analysisAlphaVY773A\'; %directory where to save input and output
% saveResults.filename = 'tracksTest1DetectionAll1.mat'; %name of file where input and output are saved
saveResults = 0; %don't save results

%verbose state
verbose = 1;

%problem dimension
probDim = 3;

%% tracking function call

% [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo(1:300),...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

% for i = 1 : 12
% for i = 1
%     movieInfoTmp((i-1)*1200+1:i*1200) = movieInfo((i-1)*1200+1:i*1200);
%     saveResults.filename = ['tracks1All_' sprintf('%02i',i) '.mat'];
%     [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfoTmp,...
%         costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
%     clear movieInfoTmp
% end

% i=6;
% movieInfoTmp((i-1)*1200+1:6800) = movieInfo((i-1)*1200+1:6800);
% saveResults.filename = ['tracks1All_' sprintf('%02i',i) '.mat'];
% [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfoTmp,...
%     costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
% clear movieInfoTmp


% for startFrame = 1 : 400 : 48000
%     endFrame = startFrame + 399;
%     saveResults.filename = ['tracks2Detection1_Frames' sprintf('%05i',startFrame) 'to' sprintf('%05i',endFrame) '.mat'];
%     disp(startFrame)
%     [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(...
%         movieInfo(startFrame:endFrame),costMatrices,gapCloseParam,...
%         kalmanFunctions,probDim,saveResults,verbose);
% end

% %%%%%%
% nFrames = 100;
% nSpotsPerFrame = 10;
% randMovementMultiplier = .01;
% driftSpeed = .05;
% angularVariation = .01;
% maxLoop = 100;
% %%%%%%
% 
% frames = 1:nFrames;
% movieInfo(frames) = struct('xCoord',[],'yCoord',[],'amp',[]);
% 
% angles = rand(nSpotsPerFrame,1) * 2 * pi;
% initCoords = rand(nSpotsPerFrame,2);
% idx = 0;
% disMat = createDistanceMatrix(initCoords,initCoords);
% disMat(disMat==0) = Inf;
% while any(disMat(:) < randMovementMultiplier) && idx < maxLoop 
%     [i,j] = find(disMat < randMovementMultiplier,1);
%     disMat = sort(disMat,2);
%     disMat = disMat(:,1:end-1);
%     disMat = mean(disMat,2);
%     if disMat(i) < disMat(j)
%         initCoords(i,:) = rand(1,2);
%     else
%         initCoords(j,:) = rand(1,2);
%     end
%     idx = idx + 1;
% end
% 
% if idx == maxLoop
%     disp('error: max loop attempts reached')
%     clear all
%     return
% end
% 
% stds = zeros(nSpotsPerFrame,1);
% amps = zeros(nSpotsPerFrame,2);
% 
% movieInfo(1).xCoord = [initCoords(:,1) stds];
% movieInfo(1).yCoord = [initCoords(:,2) stds];
% movieInfo(1).amp = amps;
% 
% for iFrame = frames(2:end)
%     angles = angles + (angularVariation * randn(nSpotsPerFrame,1));
%     movement = (randMovementMultiplier * randn(nSpotsPerFrame,2));
%     movement = movement + (driftSpeed * [cos(angles) sin(angles)]);
%     movieInfo(iFrame).xCoord = movieInfo(iFrame-1).xCoord;
%     movieInfo(iFrame).xCoord(:,1) = movieInfo(iFrame).xCoord(:,1) + movement(:,1);
%     movieInfo(iFrame).yCoord = movieInfo(iFrame-1).yCoord;
%     movieInfo(iFrame).yCoord(:,1) = movieInfo(iFrame).yCoord(:,1) + movement(:,2);
%     movieInfo(iFrame).amp = amps;
% end

% movieInfo.xCoord = [-2 0; -2 0];
% movieInfo.yCoord = [2 0; -2 0];
% movieInfo.amp = [1 0; 1 0];
% 
% movieInfo(2).xCoord = [-1 0; -1 0];
% movieInfo(2).yCoord = [1 0; -1 0];
% movieInfo(2).amp = [1 0; 1 0];
% 
% movieInfo(3).xCoord = [0 0];
% movieInfo(3).yCoord = [0 0];
% movieInfo(3).amp = [1 0];
% 
% movieInfo(4).xCoord = [1 0];
% movieInfo(4).yCoord = [0 0];
% movieInfo(4).amp = [1 0];
% 
% movieInfo(5).xCoord = [2 0; 2 0];
% movieInfo(5).yCoord = [1 0; -1 0];
% movieInfo(5).amp = [1 0; 1 0];
% 
% movieInfo(6).xCoord = [3 0; 3 0];
% movieInfo(6).yCoord = [2 0; -2 0];
% movieInfo(6).amp = [1 0; 1 0];

% movieInfo.xCoord = [-2 0; -2 0];
% movieInfo.yCoord = [2 0; -2 0];
% movieInfo.amp = [1 0; 1 0];
% 
% movieInfo(2).xCoord = [-1 0; -1 0];
% movieInfo(2).yCoord = [1 0; -1 0];
% movieInfo(2).amp = [1 0; 1 0];
% 
% movieInfo(3).xCoord = [0.01 0];
% movieInfo(3).yCoord = [0.01 0];
% movieInfo(3).amp = [1 0];
% 
% movieInfo(4).xCoord = [1 0; 1 0];
% movieInfo(4).yCoord = [1 0; -1 0];
% movieInfo(4).amp = [1 0; 1 0];
% 
% movieInfo(5).xCoord = [2 0; 2 0];
% movieInfo(5).yCoord = [2 0; -2 0];
% movieInfo(5).amp = [1 0; 1 0];

%%%%
window = 3;
mergedTime = 2;
%%%%

nFrames = (2*window) + mergedTime;
movieInfo(1:nFrames) = struct('xCoord',[],'yCoord',[],'amp',[],'zCoord',[]);
count = 0;
for iFrame = -window:-1
    count = count + 1;
    movieInfo(count).xCoord = [-iFrame 0; iFrame 0; iFrame+1.001 0; 4 0];
    movieInfo(count).yCoord = [iFrame 0; iFrame 0; iFrame 0; iFrame 0];
    movieInfo(count).amp = [0 0; 0 0; 0 0; 0 0];
    movieInfo(count).zCoord = movieInfo(count).xCoord;
end
for iFrame = 0:(mergedTime-2)
    count = count + 1;
    movieInfo(count).xCoord = [-0.001 0; 1 0; 4 0];
    movieInfo(count).yCoord = [(-0.001 + iFrame) 0; (-0.001 + iFrame) 0; iFrame+0.001 0];
    movieInfo(count).amp = [0 0; 0 0; 0 0];
    movieInfo(count).zCoord = movieInfo(count).xCoord;
end
count = count + 1;
movieInfo(count).xCoord = [0.001 0; 1 0; 4 0];
movieInfo(count).zCoord = [0.001 0; 1 0; 4 0];
movieInfo(count).yCoord = [(-0.001 + iFrame + 1) 0; (-0.001 + iFrame + 1) 0; iFrame+1.001 0];
movieInfo(count).amp = [0 0; 0 0; 0 0];
for iFrame = mergedTime:(window + mergedTime - 1)
    count = count + 1;
    movieInfo(count).xCoord = [-(-iFrame + mergedTime - 1) 0; -(iFrame - mergedTime + 1) 0; -(iFrame - mergedTime + 2) 0; 4 0];
    movieInfo(count).yCoord = [iFrame 0; iFrame 0; iFrame 0; iFrame 0];
    movieInfo(count).amp = [0 0; 0 0; 0 0; 0 0];
    movieInfo(count).zCoord = movieInfo(count).xCoord;
end

% for iFrame = -window:-1
%     count = count + 1;
%     movieInfo(count).xCoord = [-iFrame 0; iFrame 0];
%     movieInfo(count).yCoord = [iFrame 0; iFrame 0];
%     movieInfo(count).amp = [0 0; 0 0];
% end
% for iFrame = 0:(mergedTime-2)
%     count = count + 1;
%     movieInfo(count).xCoord = [-0.001 0];
%     movieInfo(count).yCoord = [(-0.001 + iFrame) 0];
%     movieInfo(count).amp = [0 0];
% end
% count = count + 1;
% movieInfo(count).xCoord = [0.001 0];
% movieInfo(count).yCoord = [(-0.001 + iFrame + 1) 0];
% movieInfo(count).amp = [0 0];
% for iFrame = mergedTime:(window + mergedTime - 1)
%     count = count + 1;
%     movieInfo(count).xCoord = [-(-iFrame + mergedTime - 1) 0; -(iFrame - mergedTime + 1) 0];
%     movieInfo(count).yCoord = [iFrame 0; iFrame 0];
%     movieInfo(count).amp = [0 0; 0 0];
% end

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);


%% ~~~ the end ~~~