%% Script to get intensity measure hazard curves (IMHCs)
clear all; close all; clc;

%% Inputs
% Specify output directory and filename
outputDir = '.\Output\';
outputFilename = 'hazCurves.mat';

% Load OpenSHA output
load(fullfile(outputDir, 'OpenSHAout.mat'));
rupRates = siteSeismicity.rupList(:,3);
nRup = siteSeismicity.numScenarios;

% Load GMPM output
load(fullfile(outputDir, 'GMPMout.mat'));
nT = length(T_all);

% Inputs for hazard curves
numIMpts = 1e3;
IMtestPts = logspace(-5,1,numIMpts)'; % Test points or thresholds for IMT when computing hazard curves

%% Compute hazard curves for H component of GM
IMHCs_geo = zeros(numIMpts,nT);
for ii=1:nT
    % Get GMPM output for current vibration period
    MUs_curr = MUs_geo_all(:,ii);
    SIGs_curr = SIGs_geo_all(:,ii);
    
    % Probability of exceedance for each rupture scenario and test point
    probExcGivenRup = 1 - logncdf(...
        repmat(IMtestPts',nRup,1),...
        repmat(MUs_curr,1,numIMpts),...
        repmat(SIGs_curr,1,numIMpts)); % nRup x numIMpts

    % Get hazard curve for current vibration period
    IMHCs_geo(:,ii) = (rupRates' * probExcGivenRup)';
end

%% Compute hazard curves for V component of GM
IMHCs_V = zeros(numIMpts,nT);
for ii=1:nT
    % Get GMPM output for current vibration period
    MUs_curr = MUs_V_all(:,ii);
    SIGs_curr = SIGs_V_all(:,ii);
    
    % Probability of exceedance for each rupture scenario and test point
    probExcGivenRup = 1 - logncdf(...
        repmat(IMtestPts',nRup,1),...
        repmat(MUs_curr,1,numIMpts),...
        repmat(SIGs_curr,1,numIMpts)); % nRup x numIMpts

    % Get hazard curve for current vibration period
    IMHCs_V(:,ii) = (rupRates' * probExcGivenRup)';
end

%% Save data
save(fullfile(outputDir,outputFilename),...
    'T_all','IMtestPts',...
    'IMHCs_geo','IMHCs_V');