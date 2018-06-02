function [fp] = dfdx(x, fx)
% -
% This file is part of KAFGA.
% Copyright (C) 2015 Markus D. Petters
%
% KAFGA is supplementary code to the following manuscript
% Petters, M. D., Kreidenweis, S. M., Ziemann, P. J.
% Prediction of cloud condensation nuclei activity for organic
% compounds using functional group contribution methods 
% Geoscientific Model Development, Manuscript number gmd-2015-172 
%
% KAFGA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% KAFGA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with KAFGA.  If not, see <http://www.gnu.org/licenses/>.
%
%- 
    nx = numel(fx);
    fp = zeros(nx,1);
    for j = 2:nx-1
        fp(j) = (fx(j+1) - fx(j-1))/(x(j+1) - x(j-1));
    end
    fp(1) = (fx(2)-fx(1))/(x(2)-x(1));
    fp(nx) = (fx(nx)-fx(nx-1))/(x(nx)-x(nx-1));
end