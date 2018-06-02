function [vw MW_C MW_H MW_O MW_N V_C V_O V_H ...
          V_N Ck Hk Ok Nk Vk Rk Qk] = load_constants()
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

sourcepath = '../../src/';  % Path to the the source code

%% funamental constants
vw = 18.07;      % molar volume of water (cm3 mol-1)
MW_C = 12.0107;  % molecular weight of carbon
MW_H = 1.00794;  % molecular weight of hydrogen
MW_O = 15.9994;  % molecular weight of oxygen
MW_N = 14.0067;  % molecular weight of nitrogen


%% Element volume parameters from Girolami
V_C = 2;         
V_O = 2;
V_H = 1;
V_N = 2;


%% load the element tables for the functional groups
Ck = load([sourcepath '/txt/Ck.txt']);
Hk = load([sourcepath '/txt/Hk.txt']);
Ok = load([sourcepath '/txt/Ok.txt']);
Nk = load([sourcepath '/txt/Nk.txt']);
Vk = load([sourcepath '/txt/Vk.txt']);

%% load the group parameters
Rk = load([sourcepath '/txt/Rk.txt']);
Qk = load([sourcepath '/txt/Qk.txt']);
end