sourcepath = '../../src/';  % Path to the the source code
addpath(sourcepath);        % Path to the the source code
cpath = 'f02comp/';         % Path to compund text files
opath1 = 'f02out/';         % Output path for .mat files 
opath2 = 'f02mat/';         % Output path for .out ASCII files
coeff = 'anm_revised.txt';  % Interaction coefficients
summary = 'summaryf02';     % Summary output file name
[s] = f02_compounds();      % Loads a structure with compound info
T = 298.15;                 % Model temperature
m = 10000; n = 1000;        % Model resolution 

myCluster = parcluster('local'); % Cluster for parallel processing
myCluster.NumWorkers = 36;       % Cluster for parallel processing
parpool(myCluster,36)            % Cluster for parallel processing

%% Main loop. Replace parfor with 'for' for regular processing             
tic
parfor i = 1:length(s)
    main(s(i), cpath, opath1, opath2, T, n, m, summary, coeff); 
end
toc