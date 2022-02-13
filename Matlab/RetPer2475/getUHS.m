%% Script to get multicomponent UHS
clear all; close all; clc;

%% Inputs
% Specify target return period
tgtRetPer = 2475; % Make sure this number agrees with current folder
rateTgt = 1/tgtRetPer; % Corresponding annual rate of exceedance

% Specify output directory and filename
outputDir = '.\Output\';
outputFilename = 'uhs.mat';

% Load hazard curves
PSHAdir = '..\PSHA\Output\';
load(fullfile(PSHAdir, 'hazCurves.mat'));

% Specify vibration periods for computing UHS
Tuhs = T_all;

%% Get UHS (horiz, RotD50)
UHS_horiz_geo = zeros(size(Tuhs));
for ii=1:length(Tuhs)
    % Set local var
    IMHCcurr = IMHCs_geo(:,ii);
    
    % ID pts closest to tgt rate
    idTgt = find( IMHCcurr<=(rateTgt*10) & IMHCcurr>=(rateTgt/10) );
    
    % Interpolate on log scale
    UHS_horiz_geo(1,ii) = exp( interp1( log(IMHCcurr(idTgt)), log(IMtestPts(idTgt)), log(rateTgt) ) );
end

%% Get UHS (vert, V data)
UHS_vert = zeros(size(Tuhs));
for ii=1:length(Tuhs)
    % Set local var
    IMHCcurr = IMHCs_V(:,ii);
    
    % ID pts closest to tgt rate
    idTgt = find( IMHCcurr<=(rateTgt*10) & IMHCcurr>=(rateTgt/10) );
    
    % Interpolate on log scale
    UHS_vert(1,ii) = exp( interp1( log(IMHCcurr(idTgt)), log(IMtestPts(idTgt)), log(rateTgt) ) );
end

%% Save data
save(fullfile(outputDir,outputFilename),...
    'tgtRetPer','rateTgt','Tuhs','UHS_horiz_geo','UHS_vert');