function [out] = unifac(c, T, n, outfile, matrix)
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

    xs = [linspace(0.00001, 0.99999, n)]; 
    out = zeros(n, 4);

    for i = 1:numel(xs)
        x = [xs(i) (1-xs(i))];
        lngR = residualActivityCoefficients(c, x, T, matrix);   
        lngC = combinatorialActivityCoefficients(c, x, T, matrix) ;
        lng = lngR + lngC;
        out(i, 1) = xs(i);
        out(i, 2) = 1-xs(i);
        out(i, 3) = exp(lng(1));
        out(i, 4) = exp(lng(2));
    end
    dlmwrite(outfile, out, 'delimiter', '\t', 'precision', 10);
end

function lngC = combinatorialActivityCoefficients(c, x, T, matrix) 
    z = 10.0;
    [R, Q, a] = load_coefficients(T, matrix);
  
    % average surface fraction and segment fraction
    r = zeros(1, numel(x));
    q = zeros(1, numel(x));

    f = @(x) x(:); % Flatten a 2D matrix
    
    for i = 1:numel(x)       
        r(i) = r(i) + sum(f(c(i).v .* R));
        q(i) = q(i) + sum(f(c(i).v .* Q));
    end

    % li
    Theta = q .* x / sum(q .* x);
    Phi = r .* x / sum(r .* x);
    l = z / 2 .* (r - q) - (r - 1);

    % activity coefficients
    lngC =  log(Phi ./ x) + z / 2 .* q .* log(Theta ./ Phi) + l - ...
        Phi ./ x * sum(x .* l);  
end

function lngR = residualActivityCoefficients(c, x, T, matrix)
    [R, Q, a] = load_coefficients(T, matrix);
    
    % Pure 1
    X = groupMoleFractions(c, [1 0]);
    Th = areaMoleFractions(c, X, Q);
    lnGk(1,:,:) = groupActivityCoefficients(c, a, Th, Q);

    % Pure 2
    X = groupMoleFractions(c, [0 1]);
    Th = areaMoleFractions(c, X, Q);
    lnGk(2,:,:) = groupActivityCoefficients(c, a, Th, Q);

    % Mixture
    X = groupMoleFractions(c, x);
    Th = areaMoleFractions(c, X, Q);
    lnGk(3,:,:) = groupActivityCoefficients(c, a, Th, Q);
   
    % For all compounds (i)
    for i = 1:numel(x)
        lngR(i) = 0;

        % For all groups (k)
        ptr = find(c(i).v > 0);
        for k = 1:numel(ptr)
            lngR(i) = lngR(i) + c(i).v(ptr(k)) * ...
                (lnGk(3, ptr(k)) - lnGk(i, ptr(k)));
        end
    end
end

function lnGk = groupActivityCoefficients(c, a, Th, Q)
    f = @(x) x(:); % Flatten a 2D matrix
    
    term1 = zeros(numel(c(1).v(:, 1)), numel(c(1).v(1, :)));
    term2 = zeros(numel(c(1).v(:, 1)), numel(c(1).v(1, :)));

    % For all groups (k)
    for n = 1:numel(c(1).v(:, 1))
        term1(n, :) = sum(f(a(:, n))' * Th);
    end

    % For all groups (k)
    for n = 1:numel(c(1).v(:, 1))
        % Evaluate sum (m)
        sum2 = 0;
        for i = 1:numel(c(1).v(:, 1))
            for j = 1:numel(c(1).v(1, :))                    
                sum2 = sum2 + a(n, i)*Th(i, j)/sum(f(a(:, i))' * Th);
            end
        end
        term2(n, :) = sum2;
    end
    
    term1;
    term2;
    lnGk = Q .* (1 - log(term1) - term2);
    lnGk(find(Th == 0)) = 0.0;    
end
        
function Th = areaMoleFractions(c, X, Q)
    Th = Q .* X / sum(Q(:) .* X(:));
end

function X = groupMoleFractions(c, x)
    f = @(x) x(:); % Flatten a 2D matrix

    % total number of groups
    sum2 = 0;
    for i = 1:numel(x)
        sum2 = sum2 + sum(f(c(i).v .* x(i)));
    end
    
    % for all groups in the mixture (preserve subgroup structure)
    for m = 1:numel(c(1).v(:,1))
        for n = 1:numel(c(1).v(1,:))
            % for all components in the mixture
            sum1 = 0; 
            for i = 1:numel(x)
                sum1 = sum1 + c(i).v(m,n)*x(i);
            end  
            X(m, n) = sum1 / sum2; 
        end
    end
end

function [Rk, Qk, a] = load_coefficients(T, matrix)
    sourcepath = '../../src/';  % Path to the the source code
    Rk = load([sourcepath '/txt/Rk.txt']);
    Qk = load([sourcepath '/txt/Qk.txt']);
    p = load(strcat([sourcepath '/txt/'], matrix));
    a = exp(-p ./ T);
end