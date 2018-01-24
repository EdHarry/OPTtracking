function visuPoints2D( coords )
%OVERLAYTRACKSONIMAGE Plot 2D Tracks
%   EHarry Nov 2014, based on plotTracks2D from uTrack-2.1.1 by KJaqaman

extendedColors = lines(32);

numTimePoints = size(coords,1);

imageInfo.coords = coords;
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
        
        if initialPlot % initially set max range
            maxX = max(imageInfo.coords(currentTimePoint,1).coords(:,1));
            minX = min(imageInfo.coords(currentTimePoint,1).coords(:,1));
            maxY = max(imageInfo.coords(currentTimePoint,1).coords(:,2));
            minY = min(imageInfo.coords(currentTimePoint,1).coords(:,2));
            subplot(1,2,1);
            set(gca,{'xlim','ylim'},{[minX maxX] , [minY maxY]})
            
            maxX = max(imageInfo.coords(currentTimePoint,2).coords(:,1));
            minX = min(imageInfo.coords(currentTimePoint,2).coords(:,1));
            maxY = max(imageInfo.coords(currentTimePoint,2).coords(:,2));
            minY = min(imageInfo.coords(currentTimePoint,2).coords(:,2));
            subplot(1,2,2);
            set(gca,{'xlim','ylim'},{[minX maxX] , [minY maxY]})
            
            subplot(1,2,1);
            hold on
            subplot(1,2,2);
            hold on
        end
        
        subplot(1,2,1);
        L = get(gca,{'xlim','ylim'});
        plot(0,0,'marker','none','YLimInclude','off','XLimInclude','off');
        set(gca,{'xlim','ylim'},L);
        
        subplot(1,2,2);
        L = get(gca,{'xlim','ylim'});
        plot(0,0,'marker','none','YLimInclude','off','XLimInclude','off');
        set(gca,{'xlim','ylim'},L);
        
        set(figH,'Name',sprintf('Frame = %2d / %2d',currentTimePoint, timeInfo.nTimePoints));
        
        %%%%%%%%%%%%%
        
        subplot(1,2,1)
        hold off
        plot(imageInfo.coords(currentTimePoint,1).coords(:,1),imageInfo.coords(currentTimePoint,1).coords(:,2),'*','LineWidth',2);
        
        subplot(1,2,2)
        hold off
        plot(imageInfo.coords(currentTimePoint,2).coords(:,1),imageInfo.coords(currentTimePoint,2).coords(:,2),'*','LineWidth',2);
       
        
        %%%%%%%%%%%%%%
    end

end

