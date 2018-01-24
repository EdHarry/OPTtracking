function params = prepPerChannelParams(params,nChan)

if isfield(params,'PerChannelParams') && ~isempty(params.PerChannelParams) && iscell(params.PerChannelParams);
      
    nPCPar = numel(params.PerChannelParams);
    
    for j = 1:nPCPar
        
        nEl = numel(params.(params.PerChannelParams{j}));
        if  nEl == 1
            params.(params.PerChannelParams{j}) = repmat(params.(params.PerChannelParams{j}),[1 nChan]);
        elseif nEl ~= nChan
            error(['The parameter "' params.PerChannelParams{j} '" was designated as a per-channel parameter, but contained ' num2str(nEl) ' elements - this must be specified as either a scalar or have array of size equal to the number of channels!'])
        end                                    
    end        
end
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
