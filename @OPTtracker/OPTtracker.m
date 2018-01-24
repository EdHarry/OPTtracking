classdef OPTtracker < tracker3D
    %OPTTRACKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        coords2D
    end
    
    methods
        function tr = OPTtracker(varargin)
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
            tracker_2D = tracker2D('coords',coords,'aor',aor,'angle',angle,'axisAngle',aAngle);
            disp('2D tracking...');
            tracker_2D.optimisedTracking;
            disp('3D object construction...');
            tracker_2D.generate3DFeatures;
            
            tr@tracker3D('coords',tracker_2D.coords3D,'aor',aor,'angle',angle,'axisAngle',aAngle);
            disp('3D tracking...');
            tr.optimisedTracking;
            tr.coords2D = coords;
            delete(tracker_2D);
        end
    end
    
end

