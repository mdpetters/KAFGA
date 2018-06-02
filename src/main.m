function [e, scun] = main(s, cpath, opath1, opath2, T, n, m, summary, ...
                           matrix)
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
t = cputime;

%% Load the lookup table for kappa and define lookup functions
load([sourcepath '/mat/lookup.mat']);
lq = @(x) min((Dlookup-x).^2.0);
kint = @(sc, in) interp1(Slookup(in,:),klookup, sc);

%% A simple terneray asignment function
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

%% Loads constansts
[vw MW_C MW_H MW_O MW_N V_C V_O V_H V_N Ck Hk Ok Nk Vk Rk Qk] = ...
    load_constants(); 

%% For all compounds
for i = 1:numel(s)

    %% Load the compound and run UNIFAC subroutine
    c(1).v = load([cpath s(i).file]);
    c(2).v = load([sourcepath '/txt/water.txt']);
    out = unifac(c, T, n, [opath1 s(i).file '.out'], matrix);

    %% Assign the activity coefficients from UNIFAC
    out = load([opath1 s(i).file '.out']);
    xw = out(:,2);      
    gs = out(:,3);
    gw = out(:,4);
    
    %% Compute the elemental composition
    f = @(x) x(:);      % Flatten a 2D matrix
    C = sum(f(Ck.*c(1).v));  % n carbon atoms
    H = sum(f(Hk.*c(1).v));  % m hydrogen atoms
    O = sum(f(Ok.*c(1).v));  % n oxygen atoms
    N = sum(f(Nk.*c(1).v));  % O nitrogen atoms
 
    %% Compute the molecular weight
    MW = C*MW_C + H*MW_H + O*MW_O + N*MW_N;  
                
    %% Density from Girolami et al. from elemental ratios
    % fc is a correction for increased densities of OH and COOH
    % max correction should be 30%. Unfused rings are not accounted for
    fc = (1+0.1*iif(sum(f(Vk.*c(1).v))<=3,sum(f(Vk.*c(1).v)),true,3));
    vs = 5*(C*V_C + H*V_H + O*V_O + N*V_N)/fc;
       
    %% Molar volume kappa with Flory Huggins correction 
    c1  = 0.014*vs^1.14 - 1.4; % Flory Huggins correction fit
    kfh = (vw + c1)/vs ;       % ideal kappa based on molar volume
    kr  = vw/vs;               % Raoult kappa based on molar volume
    
    %% Derived quantities     
    nx = numel(xw);
    aw = gw .* xw;                                % water activity
    gex = xw .* log(gw) + (1 - xw) .* log(gs);    % excess gibbs, mix
    gid = xw .* log(xw) + (1 - xw) .* log(1-xw);  % ideal gibbs, mix
    gmix = gex + gid;                             % real gibbs, mix
    dg1 = dfdx(xw, gmix);                         % derivatives 
    dg2 = dfdx(xw, dg1);
    dg3 = dfdx(xw, dg2);

    %% Phase split using area method of Eubank et al. (1992)

    % compute the area criterion by integration 
    gmixp = gmix;
    gmixp(find(gmixp > 0)) = 0; % The area method apparently fails
                                % if gmix > 1.
    for j = 1:nx
        for k = j+1:nx           
            A(j,k) = abs((gmixp(k)+gmixp(j)) * (xw(k)-xw(j))/2) ...
                - abs(trapz(xw(j:k), gmixp(j:k)));
        end
    end

    [mA,j]=max(A(:));               % find where A is max
    [j,k] = ind2sub(size(A),j);    % find where the index pair j, k

    % If Check if area is negative
    if mA <= 0
        j = numel(aw);
        k = 1;
    else 
        aw(j:k) = aw(k);         % aw is constant in the phase gap
    end
    awsat = aw(k);               % water activity of saturated sol
    xwsat = xw(j);               % composition of saturated sol
    xwsatorg = xw(k);            % composition of saturated sol org

    % sc from UNIFAC evaluated at dry diameter Dd
    Dd = s(i).D*1d-9;
    cxw = linspace(0.0001, 0.9999, m);
    caw = interp1(xw, aw, cxw);
    D = (Dd^3 .* (cxw - kr.*cxw - 1) ./ (cxw - 1)).^(1/3);
    scun = (max(caw.*exp(2.1e-9 ./ D))-1)*100;

    % find the apparent kappa at standard state from lookup table
    [~,in] = lq(Dd*1d9);
    kun = kint(scun,in);
    kun = iif(isnan(kun),1e-6,true, kun);
    
    %% write the output into the data structure
    s(i).v = c(1).v;         % group composition array
    s(i).kfh = kfh;          % molar volume kappa from Flory-Huggins
    s(i).kr = kr;            % molar volume kappa from Raoult's law
    s(i).kun = kun;          % UNIFAC kappa without conditionals 
    s(i).scun = scun;        % critical supersaturation from UNIFAC
    s(i).awsat = awsat;      % water activity of saturated solution
    s(i).xwsat = xwsat;      % mole fraction of water of sat. solution
    s(i).xwsatorg = xwsatorg;% mole fraction of water of sat. solution
    s(i).xw = xw;            % mole fraction array of water 
    s(i).gs = gs;            % activity coefficient array of organic
    s(i).gw = gw;            % activity coefficient array of water
    s(i).aw = aw;            % water activity array
    s(i).gex = gex;          % excess gibbs energy of mixing array 
    s(i).gid = gid;          % ideal gibbs energy of mixing array
    s(i).gmix = gmix;        % real gibbs energy of mixing array
    s(i).dg1 = dg1;          % first derivative dgmix/dxw
    s(i).dg2 = dg2;          % second derivative dgmix^2/d^2xw
    s(i).dg3 = dg3;          % third derivative dgmix^3/d^3xw
    s(i).vs = vs;            % molar volume of organic compound
    s(i).MW = MW;            % molecular weight of organic compound
    s(i).C = C;              % number of carbon atoms
    s(i).O = O;              % number of oxygen atoms
    s(i).H = H;              % number of hydrogen atoms
    s(i).N = N;              % number of nitrogen atoms

    %% save the individual output into a .mat file
    r = s(i);                
    file = [opath2 s(i).file '.mat'];
    save(file, 'r');

    if mod(i, 10) == 0
        i
    end
end

%% Write summary text file with important parameters
fptr = fopen(summary, 'a');
%fprintf(fptr, ['Name                            C       H       O' ...
%               '       N         MW      vs       ' ...
%               'kfh             kr             kun  ' ...
%               '            awsat           xwsat          xwsatorg']);

%fprintf(fptr, '\n', '');

for i = 1:numel(s)
    fprintf(fptr, ['%-30s\t%i\t%i\t%i\t%i\t%6.1f\t%6.1f\t%8.2e\t' ...
                   '%8.2e\t%8.2e\t%8.2e\t%8.2e\t%8.2e\t%8.2e'...
                   '\t%8.2e\t%8.2e\t%8.2e\t%8.2e'], ...
            s(i).file, s(i).C, s(i).H, s(i).O, s(i).N, s(i).MW, ...
            s(i).vs, s(i).kfh, s(i).kr, s(i).kun, s(i).awsat, ...
            s(i).xwsat, s(i).xwsatorg);
    fprintf(fptr, '\n', '');
end
fclose(fptr);

e = cputime-t;
 
end