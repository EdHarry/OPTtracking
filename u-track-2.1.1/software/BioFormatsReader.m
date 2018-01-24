classdef  BioFormatsReader < Reader
    % Concrete implementation of MovieObject for a single movie
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
    
    properties (Transient =true)
        formatReader
    end
    
    methods
        %% Constructor
        function obj = BioFormatsReader(path, iSeries)
            bfCheckJavaPath(); % Check loci-tools.jar is in the Java path
            loci.common.DebugTools.enableLogging('OFF');
            obj.formatReader = bfGetReader(path, false);
            if nargin>1 %&& ~isempty(iSeries)
                obj.formatReader.setSeries(iSeries);
%             else
%                 iSeries = obj.getSeries();
%                 obj.formatReader.setSeries(iSeries);
            end
        end
        
        function metadataStore = getMetadataStore(obj)
            r = obj.formatReader;
            metadataStore = r.getMetadataStore();
        end
        
        function series = getSeries(obj)
            series = obj.formatReader.getSeries();
        end
        
        function sizeX = getSizeX(obj, varargin)
            sizeX = obj.getMetadataStore().getPixelsSizeX(obj.getSeries()).getValue();
        end
        
        function sizeY = getSizeY(obj, varargin)
            sizeY = obj.getMetadataStore().getPixelsSizeY(obj.getSeries()).getValue();
        end
        
        function sizeZ = getSizeZ(obj, varargin)
            sizeZ = obj.getMetadataStore().getPixelsSizeZ(obj.getSeries()).getValue();
        end
        
        function sizeT = getSizeT(obj, varargin)
            sizeT = obj.getMetadataStore().getPixelsSizeT(obj.getSeries()).getValue();
        end
        
        function sizeC = getSizeC(obj, varargin)
            sizeC = obj.getMetadataStore().getPixelsSizeC(obj.getSeries()).getValue();
        end
        
        function bitDepth = getBitDepth(obj, varargin)
            pixelType = obj.formatReader.getPixelType();
            bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
            bitDepth = 8 * bpp;
        end
        
        function fileNames = getImageFileNames(obj, iChan, varargin)
            % Generate image file names
            [~, fileName] = fileparts(char(obj.formatReader.getCurrentFile));
            basename = sprintf('%s_s%g_c%d_t',fileName, obj.getSeries()+1, iChan);
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
        end
        
        function channelNames = getChannelNames(obj, iChan)
            [~, fileName, fileExt] = fileparts(char(obj.formatReader.getCurrentFile));
            
            if obj.formatReader.getSeriesCount() > 1
                base = [fileName fileExt ' Series ' num2str(obj.getSeries()+1) ' Channel '];
            else
                base = [fileName fileExt ' Channel '];
            end
            
            channelNames = arrayfun(@(x) [base num2str(x)], iChan, 'Unif',false);
        end
        
        function index = getIndex(obj, z, c, t)
            index = loci.formats.FormatTools.getIndex(obj.formatReader, z, c, t);
        end
        
        function I = loadImage(obj, c, t)
            
            ip = inputParser;
            ip.addRequired('c', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            ip.addRequired('t', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT()));
            ip.parse(c, t);
            
            % Using bioformat tools, get the reader and retrieve dimension order
            I = bfGetPlane(obj.formatReader, obj.getIndex(0, c-1, t-1) + 1);
        end
        
        function delete(obj)
            obj.formatReader.close()
        end
    end
end