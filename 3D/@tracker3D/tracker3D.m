classdef tracker3D < handle
    %TRACKER2D 2D tracker
    %   EHarry Nov 2014
    
    % main properties
    properties (SetAccess = private)
        axisOfRotation
        rotationAngle
        axisAngle
        rawCoords
        alignedCoords
        trackingCoords
        tracks
    end
    
    % dependent properties
    properties (SetAccess = private, Dependent = true)
        tracks_raw
        tracks_tracking
        tracks_filtered
        tracks_filtered_raw
        tracks_filtered_tracking
    end
    
    % main methods
    methods
        function tr = tracker3D(varargin) % constructor
            %% SETUP
            
            % defaults
            param = struct('coords',[],'aor',[0 0],'angle',0,'axisAngle',pi/2);
            
            % user-selected param
            var = fieldnames(param);
            for i = 1:2:length(varargin)
                if ~isempty(strcmpi(varargin{i},var))
                    param.(varargin{i}) = varargin{i+1};
                end
            end
            
            % unpack param
            coords = param.coords;
            aor = param.aor;
            angle = param.angle;
            aAngle = param.axisAngle;
            
            %% CONSTRUCTOR MAIN
            
            if isempty(coords)
                error('tracker3D:noInput','must supply coordinates to tracker3D')
            end
            
            coordsA = alignFrames_3Drotation(coords, aor, angle, aAngle);
            
            tr.axisOfRotation = aor;
            tr.rotationAngle = angle;
            tr.rawCoords = coords;
            tr.alignedCoords = coordsA;
            tr.axisAngle = aAngle;
        end
        
        function optimisedTracking(tr)
            [tC, t] = tr.optimise3Dtracking(tr.alignedCoords);
            
            if ~isempty(t)
                tr.trackingCoords = tC;
                tr.tracks = t;
            end
        end
        
        function visualise(tr,selectString)
            if isempty(tr.tracks)
                disp('no tracking has been performed yet');
                return
            end
            
            if nargin < 2
                selectString = [];
            end
            
            tr.visualiser3D(selectString);
        end
        
    end
    
    % get methods
    methods
        function value = get.tracks(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.tracks, tr.alignedCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks_raw(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.tracks, tr.rawCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks_tracking(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.tracks, tr.trackingCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks_filtered(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.filter3DTracks(tr.tracks_tracking), tr.alignedCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks_filtered_raw(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.filter3DTracks(tr.tracks_tracking), tr.rawCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks_filtered_tracking(tr)
            if ~isempty(tr.tracks)
                value = tr.generateTrackingCoords(tr.filter3DTracks(tr.tracks_tracking), tr.trackingCoords);
            else
                value = [];
            end
        end
    end
    
    % hidden static methods
    methods (Hidden = true, Static = true)
        [argout1, argout2] = optimise3Dtracking(argin);
        argout = generateTrackingCoords(argin1, argin2);
        argout = filter3DTracks(argin);
    end
    
    % hidden methods
    methods (Hidden = true)
        visualiser3D(argin1, argin2);
    end
    
end

