sourcepath = '../../src/';  % Path to the the source code
addpath(sourcepath);        % Path to the the source code
cpath = 'fs01comp/';        % Path to compund text files
opath1 = 'fs01out/';        % Output path for .mat files 
opath2 = 'fs01mat/';        % Output path for .out ASCII files
coeff = 'anm_standard.txt'; % Interaction coefficients
summary = 'summaryfs01';    % Summary output file name
[s] = fs01_compounds();     % Loads a structure with compound info
T = 298.15;                 % Model temperature
m = 10000; n = 1000;        % Model resolution 

tic
for i = 1:length(s)
    main(s(i), cpath, opath1, opath2, T, n, m, summary, coeff); 
end
toc

%% Generates the supplementary figure                                  
set(gcf, 'PaperUnits', 'inches','PaperSize', [6.6 3]); 
set(gcf, 'PaperPosition', [0 0 6.6 3], 'Color', 'w');

load('fs01mat/ethanol.txt.mat')
l = read_ames('/ethanol_water.txt');

subplot(121)
plot(1-r.xw, r.gw, 'k')
hold on
plot(l.x1, l.g2, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4)
plot(1-r.xw, r.gs, 'r')
plot(l.x1, l.g1, 'ro', 'MarkerFaceColor', 'w', 'MarkerSize', 4)
hold off

xlabel('Mole fraction (-)'); ylabel('Activity coefficient (-)');
set(gca, 'FontSize', 10, 'TickLength', [0.02 0.02], 'Box', 'on');         
set(gca, 'XminorTick', 'on', 'YminorTick', 'on');     
xlim([-0.02 1.02]); ylim([0.8, 8]);

load('fs01mat/acetone.txt.mat')
l = read_ames('/acetone_water.txt');
subplot(122)
plot(1-r.xw, r.gw, 'k')
hold on
plot(l.x1, l.g2, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4)
plot(1-r.xw, r.gs, 'b')
plot(l.x1, l.g1, 'bo','MarkerFaceColor', 'w', 'MarkerSize', 4)
hold off
xlabel('Mole fraction (-)'); ylabel('Activity coefficient (-)');
set(gca, 'FontSize', 10, 'TickLength', [0.02 0.02], 'Box', 'on');         
set(gca, 'XminorTick', 'on', 'YminorTick', 'on');     
xlim([-0.02 1.02]);ylim([0.8, 12]);

print('-dpdf','fs01.pdf')