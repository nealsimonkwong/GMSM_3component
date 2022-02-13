%% Script to save output from OpenSHA command line application into Matlab
%{
- Conduct PSHA via OpenSHA command line app before using this script
- Extract OpenSHA output into Matlab to facilitate GMSM
- By Neal (Simon) Kwong; nealsimonkwong@berkeley.edu
%}
clear all; close all; clc;

%% Inputs
% Specify site info 
siteLabel = 'DavisCA';
siteLatLong = [38.54 -121.74];

%{
NEHRP 2003 site classifications (m/sec)
A   Vs30 > 1500
B   760 < Vs30 <= 1500
C   360 < Vs30 <= 760
D   180 <= Vs30 <= 360
E   Vs30 < 180 (or with combo of low blow counts and undrained shear
strength, or high plasticity index)
F   Soils susceptible to liquefaction
%}
siteVS30 = 760; % m/s
siteNEHRPclass = 'C';

%{
Default values in command line app as follows. To maintain consistency btw 
command line and GUI apps in OpenSHA, default values in command line app
used.
%}
Z_2pt5 = 1.0; % km
Z_1pt0 = 0.100; % km

% GMPM of interest
GMPMstr = 'CB2014';

% Vibration periods of interest
T = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0];
nPer = length(T);

% Where output from OpenSHA command line app is stored
subFolderInOpenSHAdir = siteLabel;
OpenSHAdir = fullfile('..\..\OpenSHA\',subFolderInOpenSHAdir);

% Specify output directory and filename
outputDir = '.\Output';
outputFilename = ['OpenSHAout.mat'];

%% Extract ERF from OpenSHA output into variable rupList
%{
% Format of variable rupList: nRup x 6
% Col 1: SourceID
% Col 2: RuptureID
% Col 3: Rupture rate
% Col 4: Magnitude
% Col 5: Closest dist
% Col 6: J-B dist
%}

% Get rupture metadata
fileID = fopen(fullfile(OpenSHAdir,'src_rup_metadata.txt'));
C = textscan(fileID, '%f %f %f %f %*[^\n]'); % SrcID RupIDgivenSrc RupRate Mag
fclose(fileID);
rupList = [C{1} C{2} C{3} C{4}];
rupRates = C{3};
nRup = length(rupRates); % Number of rupture scenarios

% Extract closest distances
fileID = fopen(fullfile(OpenSHAdir,'rup_dist_info.txt'));
C = textscan(fileID, '%f %f %f');
fclose(fileID);
rupList(:,5) = C{3};

% Extract Joyner-Boore distances
fileID = fopen(fullfile(OpenSHAdir,'rup_dist_jb_info.txt'));
C = textscan(fileID, '%f %f %f');
fclose(fileID);
rupList(:,6) = C{3};

%% Extract GMPM output from OpenSHA
% Create strings for vibration periods
Tstr = cellstr( num2str(T',2) );
for ii=1:nPer % Revise strings for periods 1 to 10 sec and remove leading white space
    if T(ii)>=1
        Tstr{ii} = num2str(T(ii),'%2.1f');
    end
    Tstr{ii} = strrep(Tstr{ii},' ','');
end

% Get data
MUs_openSHA = zeros(nRup,nPer); SIGs_openSHA = zeros(size(MUs_openSHA));
for ii = 1:nPer
    Tcurr = T(ii);
    fnCurr = [GMPMstr '_SA_' Tstr{ii} '.txt'];
    dataCurr = load(fullfile(OpenSHAdir,fnCurr),'r');
    MUs_openSHA(:,ii) = dataCurr(:,3);
    SIGs_openSHA(:,ii) = dataCurr(:,4);
end

%% Save data into MATLAB
% Store site data
siteInfo.siteName = siteLabel;
siteInfo.lat_long = siteLatLong; % 1x2 row vector of coordinates
siteInfo.vs30 = siteVS30;
siteInfo.z2pt5 = Z_2pt5;
siteInfo.z1pt0 = Z_1pt0;

% Store ERF data
siteSeismicity.rupList = rupList; % Matrix of rupture scenarios
siteSeismicity.numEQsources = length(unique(rupList(:,1))); % Number of EQ sources
siteSeismicity.numScenarios = size(rupList,1); % Each rupture scenario is a triplet of EQ source, M, and R
siteSeismicity.EQoccRate = sum(rupList(:,3)); % Annual rate of any rupture occurring

% Save
save(fullfile(outputDir,outputFilename),'siteInfo','siteSeismicity',...
    'GMPMstr','T','MUs_openSHA','SIGs_openSHA');