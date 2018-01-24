function visuTracks2D( trackedFeatureInfo, useMatlabVersion )
%OVERLAYTRACKSONIMAGE Plot 2D Tracks
%   EHarry Nov 2014, based on plotTracks2D from uTrack-2.1.1 by KJaqaman

if nargin < 2
    useMatlabVersion = false;
end

extendedColors = lines(32);

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    if isfield(trackedFeatureInfo, 'coords')
        trackedFeatureInfo = constructTrackMatFromCoords(trackedFeatureInfo);
        [numTracks,numTimePoints] = size(trackedFeatureInfo);
        numTimePoints = numTimePoints/8;
    else
        tmp = vertcat(trackedFeatureInfo.seqOfEvents);
        numTimePoints = max(tmp(:,1));
        clear tmp
    end
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
end

if isstruct(trackedFeatureInfo) %if tracks are input in structure format
    
    % use Icy if selected
    if ~useMatlabVersion
        icyInterface = IcyInterface;
        icyInterface.startIcy;
        icyInterface.displayTracks(convStruct2MatIgnoreMS(trackedFeatureInfo));
        return
    end
    
    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end
    
    %if all tracks have only one segment ...
    if max(numSegments) == 1
        
        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = false;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';
        
        %store tracks in a matrix
        trackedFeatureInfo = NaN(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end
        
    else %if some tracks have merging/splitting branches
        
        %indicate that in the variable mergeSplit
        mergeSplit = true;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
        end
        
        %put all tracks together in a matrix
        trackedFeatureInfo = NaN(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end
        
    end
    
else %if tracks are not input in structure format
    
    % use Icy if selected
    if ~useMatlabVersion
        icyInterface = IcyInterface;
        icyInterface.startIcy;
        icyInterface.displayTracks(trackedFeatureInfo);
        return
    end
    
    inputStructure = [];
    
    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = false;
    
    %indicate that each track consists of one segment
    numSegments = ones(numTracks,1);
    
    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';
    
end

tracksX = trackedFeatureInfo(:,1:8:end)';
tracksY = trackedFeatureInfo(:,2:8:end)';

imageInfo.inputStructure = inputStructure;
imageInfo.numTracks = numTracks;
imageInfo.tracksX = tracksX;
imageInfo.tracksY = tracksY;
imageInfo.mergeSplit = mergeSplit;
imageInfo.trackStartRow = trackStartRow;
imageInfo.numSegments = numSegments;
imageInfo.extendedColors = extendedColors;

% initialize plot window
figH = figure('Name',sprintf('Frame = %2d / %2d',1, numTimePoints));
set(figH,'KeyPressFcn',@figure_keyPress);
% user data of the figure stores the current time point and the number of frames
imageInfo.timeInfo = struct('nTimePoints',numTimePoints,...
    'currentTimePoint',1);
set(figH,'UserData',imageInfo);
plotTracksLocal(figH,true);


%% LOCAL FUNCTIONS

    function figure_keyPress(src,event)
        
        imageInfo = get(src,'UserData');
        timeInfo = imageInfo.timeInfo;
        
        switch event.Key
            case 'rightarrow'
                timeInfo.currentTimePoint = timeInfo.currentTimePoint + 1;
            case 'leftarrow'
                timeInfo.currentTimePoint = timeInfo.currentTimePoint - 1;
        end
        
        if timeInfo.currentTimePoint > timeInfo.nTimePoints
            timeInfo.currentTimePoint = 1;
        end
        if timeInfo.currentTimePoint < 1
            timeInfo.currentTimePoint = timeInfo.nTimePoints;
        end
        
        imageInfo.timeInfo = timeInfo;
        
        set(src,'UserData',imageInfo);
        plotTracksLocal(src,false);
        
    end % end of function figure_keyPress

    function plotTracksLocal(figH,initialPlot)
        
        % figure is still open
        figure(figH)
        % make invisibe while processing
        %         set(figH,'Visible','off');
        %cameraProps = get(gca,cameraPropertyNames);
        
        imageInfo = get(figH,'UserData');
        timeInfo = imageInfo.timeInfo;
        currentTimePoint = timeInfo.currentTimePoint;
        timeRange = [max([1 (currentTimePoint-10)]) currentTimePoint];
        
        if initialPlot % initially set max range
            maxX = max(imageInfo.tracksX(:));
            minX = min(imageInfo.tracksX(:));
            maxY = max(imageInfo.tracksY(:));
            minY = min(imageInfo.tracksY(:));
            set(gca,{'xlim','ylim'},{[minX maxX] , [minY maxY]})
            hold on
        end
        
        L = get(gca,{'xlim','ylim'});
        plot(0,0,'marker','none','YLimInclude','off','XLimInclude','off');
        set(figH,'Name',sprintf('Frame = %2d / %2d',currentTimePoint, timeInfo.nTimePoints));
        set(gca,{'xlim','ylim'},L)
        %         if ~initialPlot
        %             zoom reset
        %             set(gca,{'xlim','ylim'},L)
        %         end
        
        %%%%%%%%%%%%%%
        
        %calculate the number of time points to be plotted
        numTimePlot = timeRange(2) - timeRange(1) + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plotting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        hold on
        
        %extract the portion of tracksX and tracksY that is of interest
        tracksXP = imageInfo.tracksX(timeRange(1):timeRange(2),:);
        tracksYP = imageInfo.tracksY(timeRange(1):timeRange(2),:);
        
        colorTime = '3';
        markerType = 'none';
        
        switch colorTime
            
            case '1' %if user wants to color-code time
                
                xData=arrayfun(@(x)[tracksXP(~isnan(tracksXP(:,x)),x); NaN],1:size(tracksXP,2),'Unif',false);
                yData=arrayfun(@(x)[tracksYP(~isnan(tracksYP(:,x)),x); NaN],1:size(tracksYP,2),'Unif',false);
                plot(vertcat(xData{:}),vertcat(yData{:}),'k:','YLimInclude','off','XLimInclude','off');
                
                %get the overall color per time interval
                colorOverTime = timeColormap(numTimePlot);
                
                %overlay tracks with color coding wherever a feature has been detected
                for i=1:numTimePlot-1
                    validData=~all(isnan(tracksXP(i:i+1,:)),1);
                    xData=vertcat(tracksXP(i:i+1,validData),NaN(1,sum(validData)));
                    yData=vertcat(tracksYP(i:i+1,validData),NaN(1,sum(validData)));
                    plot(xData(:),yData(:),'color',colorOverTime(i,:),'YLimInclude','off','XLimInclude','off');
                end
                
            case '2' %no time color-coding, loop through series of colors to color tracks
                
                %plot tracks by looping through colors
                %missing intervals are indicated by a dotted line
                for i = 1 : trackStartRow(end) + numSegments(end) - 1
                    obsAvail = find(~isnan(tracksXP(:,i)));
                    plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:','YLimInclude','off','XLimInclude','off');
                    plot(tracksXP(:,i),tracksYP(:,i),'color',colorLoop(mod(i-1,7)+1,:),...
                        'marker',markerType,'YLimInclude','off','XLimInclude','off');
                end
                
            case '3' % no time color-coding, use extendedColors
                
                %plot tracks by looping through colors
                %missing intervals are indicated by a dotted line
                for i = 1 : trackStartRow(end) + numSegments(end) - 1
                    obsAvail = find(~isnan(tracksXP(:,i)));
                    plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:','YLimInclude','off','XLimInclude','off');
                    plot(tracksXP(:,i),tracksYP(:,i),'color',imageInfo.extendedColors(mod(i-1,23)+1,:),...
                        'marker',markerType,'YLimInclude','off','XLimInclude','off');
                    % put marker at current position
                    plot(tracksXP(end,i),tracksYP(end,i),'color',imageInfo.extendedColors(mod(i-1,23)+1,:),...
                        'marker','+','YLimInclude','off','XLimInclude','off');
                end
                
            otherwise %no time color-coding, all tracks same color
                
                %plot tracks with the line color indicated
                %missing intervals are indicated by a dotted line
                for i = 1 : imageInfo.trackStartRow(end) + imageInfo.numSegments(end) - 1
                    obsAvail = find(~isnan(tracksXP(:,i)));
                    plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:','YLimInclude','off','XLimInclude','off');
                    plot(tracksXP(:,i),tracksYP(:,i),colorTime,'marker',markerType,'YLimInclude','off','XLimInclude','off');
                    %             plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),':','Color',[0.7 0.7 0.7]);
                    %             plot(tracksXP(:,i),tracksYP(:,i),'marker',markerType,'Color',[0.7 0.7 0.7]);
                end
                
        end %(switch colorTime)
        
        %show merges and splits
        if imageInfo.mergeSplit
            
            %go over all tracks
            for iTrack = 1 : imageInfo.numTracks
                
                %parse sequence of events of this compound track and find merges and
                %splits
                seqOfEvents = imageInfo.inputStructure(iTrack).seqOfEvents;
                indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) ...
                    & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
                indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) ...
                    & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
                
                %go over all splits
                for iSplit = indxSplit
                    
                    %get time of splitting
                    timeSplit = seqOfEvents(iSplit,1);
                    
                    %determine row where starting track is located
                    rowS = imageInfo.trackStartRow(iTrack) + seqOfEvents(iSplit,3) - 1;
                    
                    %determine row where splitting track is located
                    rowSp = imageInfo.trackStartRow(iTrack) + seqOfEvents(iSplit,4) - 1;
                    
                    %plot split as a dash-dotted line
                    plot([imageInfo.tracksX(timeSplit,rowS) imageInfo.tracksX(timeSplit-1,rowSp)], ...
                        [imageInfo.tracksY(timeSplit,rowS) imageInfo.tracksY(timeSplit-1,rowSp)],'k-.','YLimInclude','off','XLimInclude','off');
                    
                end
                
                %go over all merges
                for iMerge = indxMerge
                    
                    %get time of merging
                    timeMerge = seqOfEvents(iMerge,1);
                    
                    %determine row where ending track is located
                    rowE = imageInfo.trackStartRow(iTrack) + seqOfEvents(iMerge,3) - 1;
                    
                    %determine row where merging track is located
                    rowM = imageInfo.trackStartRow(iTrack) + seqOfEvents(iMerge,4) - 1;
                    
                    %plot merge as a dashed line
                    plot([imageInfo.tracksX(timeMerge-1,rowE) imageInfo.tracksX(timeMerge,rowM)], ...
                        [imageInfo.tracksY(timeMerge-1,rowE) imageInfo.tracksY(timeMerge,rowM)],'k--','YLimInclude','off','XLimInclude','off');
                    %             plot([tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                    %                 [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'--','Color',[0.7 0.7 0.7]);
                    
                end
                
            end %(for iTrack = 1 : numTracks)
            
        end %(if mergeSplit)
        
        indicateSE = true;
        
        if indicateSE %if user wants to indicate starts and ends
            
            %if there are merges and splits
            if imageInfo.mergeSplit
                
                %go over all tracks
                for iTrack = 1 : imageInfo.numTracks
                    
                    %parse sequence of events of this compound track and find starts and
                    %ends
                    seqOfEvents = imageInfo.inputStructure(iTrack).seqOfEvents;
                    indxStart = (find(seqOfEvents(:,2) == 1 & isnan(seqOfEvents(:,4)) ...
                        & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
                    indxEnd = (find(seqOfEvents(:,2) == 2 & isnan(seqOfEvents(:,4)) ...
                        & seqOfEvents(:,1) >= timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
                    
                    %get the information of the starts
                    startInfo = [];
                    for i = 1 : length(indxStart)
                        iStart = indxStart(i);
                        
                        %get start time
                        timeStart = seqOfEvents(iStart,1);
                        
                        %determine row where starting track is located in big matrix
                        %of tracks
                        rowS = imageInfo.trackStartRow(iTrack) + seqOfEvents(iStart,3) - 1;
                        
                        %get coordinates at the start
                        startInfo(i,:) = [imageInfo.tracksX(timeStart,rowS) imageInfo.tracksY(timeStart,rowS) timeStart];
                        
                    end
                    
                    %get the information of the ends
                    endInfo = [];
                    for i = 1 : length(indxEnd)
                        iEnd = indxEnd(i);
                        
                        %get end time
                        timeEnd = seqOfEvents(iEnd,1);
                        
                        %determine row where ending track is located in big matrix
                        %of tracks
                        rowE = imageInfo.trackStartRow(iTrack) + seqOfEvents(iEnd,3) - 1;
                        
                        %get coordinates at the end
                        endInfo(i,:) = [imageInfo.tracksX(timeEnd,rowE) imageInfo.tracksY(timeEnd,rowE) timeEnd];
                        
                    end
                    
                    %place circles at track starts and squares at track ends
                    switch colorTime
                        case {'1','2','3'}
                            if ~isempty(startInfo)
                                plot(startInfo(:,1),startInfo(:,2),'k',...
                                    'LineStyle','none','marker','o','YLimInclude','off','XLimInclude','off');
                            end
                            if ~isempty(endInfo)
                                plot(endInfo(:,1),endInfo(:,2),'k',...
                                    'LineStyle','none','marker','square','YLimInclude','off','XLimInclude','off');
                            end
                            % JD 2/09: case 2 is identical to case 1
                            %                 case '2'
                            %                     if ~isempty(startInfo)
                            %                         plot(startInfo(:,1),startInfo(:,2),'k',...
                            %                             'LineStyle','none','marker','o');
                            %                     end
                            %                     if ~isempty(endInfo)
                            %                         plot(endInfo(:,1),endInfo(:,2),'k',...
                            %                             'LineStyle','none','marker','square');
                            %                     end
                        otherwise
                            if ~isempty(startInfo)
                                plot(startInfo(:,1),startInfo(:,2),colorTime,...
                                    'LineStyle','none','marker','o','YLimInclude','off','XLimInclude','off');
                            end
                            if ~isempty(endInfo)
                                plot(endInfo(:,1),endInfo(:,2),colorTime,...
                                    'LineStyle','none','marker','square','YLimInclude','off','XLimInclude','off');
                            end
                    end
                    
                end %(for iTrack = 1 : numTracks)
                
            else %if there are no merges and splits
                
                %find the beginning and end of each track
                for i=imageInfo.numTracks:-1:1
                    timePoint = find(~isnan(imageInfo.tracksX(:,i)));
                    startInfo(i,:) = [imageInfo.tracksX(timePoint(1),i) ...
                        imageInfo.tracksY(timePoint(1),i) timePoint(1)];
                    endInfo(i,:) = [imageInfo.tracksX(timePoint(end),i) ...
                        imageInfo.tracksY(timePoint(end),i) timePoint(end)];
                end
                
                %place circles at track starts and squares at track ends if they happen to
                %be in the plotting region of interest
                switch colorTime
                    case {'1','2','3'}
                        indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                        plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o','YLimInclude','off','XLimInclude','off');
                        indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                        plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square','YLimInclude','off','XLimInclude','off');
                        %             case '2'
                        %                 indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                        %                 plot(startInfo(indx,1),startInfo(indx,2),'k','LineStyle','none','marker','o');
                        %                 indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                        %                 plot(endInfo(indx,1),endInfo(indx,2),'k','LineStyle','none','marker','square');
                    otherwise
                        indx = find(startInfo(:,3)>=timeRange(1) & startInfo(:,3)<=timeRange(2));
                        plot(startInfo(indx,1),startInfo(indx,2),colorTime,...
                            'LineStyle','none','marker','o','YLimInclude','off','XLimInclude','off');
                        indx = find(endInfo(:,3)>=timeRange(1) & endInfo(:,3)<=timeRange(2));
                        plot(endInfo(indx,1),endInfo(indx,2),colorTime,...
                            'LineStyle','none','marker','square','YLimInclude','off','XLimInclude','off');
                end
                
            end %(if mergeSplit)
            
        end %(if indicateSE)
        
        
        hold off
        
        %%%%%%%%%%%%%%
    end

end

