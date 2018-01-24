classdef tracker2D < handle
    %TRACKER2D 2D tracker
    %   EHarry Nov 2014
    
    % main properties
    properties (SetAccess = private)
        aor
        angle
        axisAngle
        cellCoords
        alignedCoords
        tracks2D
        tracks3D
        coords3D
    end
    
    % dependent properties
    properties (SetAccess = private, Dependent = true)
        tracks2D_aligned
        split_tracks2D
        split_tracks2D_aligned
    end
    
    % hidden properties
    properties (SetAccess = private, GetAccess = private, Hidden = true)
        tracks3D_2Dindexes
    end
    
    % main methods
    methods
        function tr = tracker2D(varargin) % constructor
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
            aor_ = param.aor;
            angle_ = param.angle;
            axisAngle_ = param.axisAngle;
            
            %% CONSTRUCTOR MAIN
            
            if isempty(coords)
                error('tracker2D:noInput','must supply coordinates to tracker2D')
            end
            
            tr.cellCoords = coords;
            tr.aor = aor_;
            tr.angle = angle_;
            tr.axisAngle = axisAngle_;
        end
        
        function optimisedTracking(tr)
            [aC, t2D] = tr.optimise2Dtracking(tr.cellCoords);
            
            if ~isempty(t2D)
                tr.alignedCoords = aC;
                tr.tracks2D = t2D;
            end
        end
        
        function visualise(tr,selectString)
            if isempty(tr.tracks2D)
                disp('no tracking has been performed yet');
                return
            end
            
            if nargin < 2
                selectString = [];
            end
            
            tr.visualiser2D(selectString);
        end
        
        function generate3DFeatures(tr)
            if isempty(tr.tracks2D)
                disp('no tracking has been performed yet');
                return
            end
            
            [tracks3D_,coords3D_,trackMatFull] = tr.generate3DTracks(tr.split_tracks2D,tr.cellCoords);
            tracks3D_ = tr.generateTrackingCoords(tracks3D_, coords3D_);
            [tr.tracks3D,tr.coords3D,tr.tracks3D_2Dindexes] = tr.augment3DCoords(tracks3D_,coords3D_,trackMatFull,tr.cellCoords,tr.aor,tr.angle,tr.axisAngle);
        end
        
    end
    
    % get methods
    methods
        function value = get.tracks2D(tr)
            if ~isempty(tr.tracks2D)
                value = tr.generateTrackingCoords(tr.tracks2D, tr.cellCoords);
            else
                value = [];
            end
        end
        
        function value = get.tracks2D_aligned(tr)
            if ~isempty(tr.tracks2D)
                value = tr.generateTrackingCoords(tr.tracks2D, tr.alignedCoords);
            else
                value = [];
            end
        end
        
        function value = get.split_tracks2D(tr)
            if ~isempty(tr.tracks2D)
                value = tr.generateTrackingCoords(tr.trackSpliter_2D(tr.tracks2D), tr.cellCoords);
            else
                value = [];
            end
        end
        
        function value = get.split_tracks2D_aligned(tr)
            if ~isempty(tr.tracks2D)
                value = tr.generateTrackingCoords(tr.trackSpliter_2D(tr.tracks2D), tr.alignedCoords);
            else
                value = [];
            end
        end
    end
    
    % hidden static methods
    methods (Hidden = true, Static = true)
        [argout1, argout2] = optimise2Dtracking(argin);
        argout = generateTrackingCoords(argin1, argin2);
        argout = trackSpliter_2D(argin);
        [argout1, argout2, argout3] = generate3DTracks(argin1, argin2);
        [argout1, argout2, argout3] = augment3DCoords(argin1, argin2, argin3, argin4, argin5, argin6, argin7);
    end
    
    % hidden methods
    methods (Hidden = true)
        visualiser2D(argin1, argin2);
    end
    
end

