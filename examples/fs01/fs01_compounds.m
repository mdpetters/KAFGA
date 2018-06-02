function [s] = fs01_compounds()
    s(1).name    = 'ethanol';
    s(1).file    = 'ethanol.txt';
    s(1).sol     = 1;
    s(1).MW_obs  = 46.06844;
    s(1).rho_obs = 0.7890;
    s(1).Cx      = 2;
    s(1).Hx      = 6;
    s(1).Ox      = 1;
    s(1).Nx      = 0;
    s(1).D       = 100;

    s(2).name    = 'acetone';
    s(2).file    = 'acetone.txt';
    s(2).sol     = 1;
    s(2).MW_obs  = 58.08;
    s(2).rho_obs = 0.791;
    s(2).Cx      = 2;
    s(2).Hx      = 6;
    s(2).Ox      = 1;
    s(2).Nx      = 0;
    s(2).D       = 100;
end

