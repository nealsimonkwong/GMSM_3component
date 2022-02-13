%% Script to get GMPM output for all rupture scenarios
% For convenience of applying GMPMs, assume all rupture scenarios in the
% ERF are vertical strike-slip 
clear all; close all; clc;

%% Inputs
% Specify output directory and filename
outputDir = '.\Output\';
outputFilename = 'GMPMout.mat';

% Load OpenSHA output
load(fullfile(outputDir, 'OpenSHAout.mat'));
nRup = siteSeismicity.numScenarios;

% Specify labels for GMPMs adopted
GMPMstr_horiz = 'CB2014'; % Campbell and Bozorgnia 2014 (H component)
GMPMstr_vert = 'BC2016_V'; % Bozorgnia and Campbell 2016 (V component)

%% Specify vibration periods of interest for computing hazard curves and UHSs
T_CB2014 = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0]; % Periods used in GMPM development
T_all = T_CB2014; % Can include conditioning periods in T_all if needed
nT = length(T_all);

%% Assumptions for each rupture scenario in ERF
FRV = 0; FNM = 0;
dip = 90; 
lambda = 0; Fhw = 0; 
region = 0;
Sj = 0;
W = 15;
Ztor = 0;
Zhyp = 10;
Zbot = 15;
Vs30 = siteInfo.vs30;

%% Get GMPM output (horizontal, RotD50)
% Initialize
MUs_geo_all = zeros(nRup,nT);
SIGs_geo_all = zeros(nRup,nT);
% Apply GMPM for all rupture scenarios
for ii=1:nRup
    [Sa,sigma] = CB_2014_nga(...
        siteSeismicity.rupList(ii,4),... % Magnitude
        T_all,...
        siteSeismicity.rupList(ii,5),... % Rrup
        siteSeismicity.rupList(ii,5),... % Rjb
        siteSeismicity.rupList(ii,5),... % Rx
        W, Ztor, Zbot, dip, lambda, Fhw,...
        siteInfo.vs30,...
        siteInfo.z2pt5,...
        Zhyp, region);
    MUs_geo_all(ii,:) = log(Sa);
    SIGs_geo_all(ii,:) = sigma;
end

%% Get GMPM output (vertical)
% Initialize
MUs_V_all = zeros(nRup,nT);
SIGs_V_all = zeros(nRup,nT);
% Apply GMPM for all rupture scenarios
for ii=1:nRup
    [Sa, sigma] = gmpmV_BC2016(...
        siteSeismicity.rupList(ii,4),... % Magnitude
        siteSeismicity.rupList(ii,5),... % Rrup
        siteSeismicity.rupList(ii,5),... % Rjb
        siteSeismicity.rupList(ii,5),... % Rx
        FRV, FNM, dip,...
        siteInfo.vs30,...
        region, Sj,...
        T_all,...
        W, Ztor,...
        siteInfo.z2pt5,...
        Zhyp);
    MUs_V_all(ii,:) = log(Sa);
    SIGs_V_all(ii,:) = sigma;    
end

%% Save GMPM output
save(fullfile(outputDir,outputFilename),...
    'GMPMstr_horiz','GMPMstr_vert',...
    'T_CB2014','T_all',...
    'MUs_geo_all','SIGs_geo_all','MUs_V_all','SIGs_V_all');