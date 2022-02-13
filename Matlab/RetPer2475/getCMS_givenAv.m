%% Get multicomponent CMS given A_V at specified T*
clear all; close all; clc;

%% Input
% Define component of GM for conditioning
compStr = 'Av';

% Specify conditioning period (for convenience, this period should have
% been used to compute hazard curves)
Tstar = 0.1; % Example periods: 0.1 0.5 0.75

% Create label for chosen T*
TstarStr = num2str(Tstar);
TstarStr = strrep(TstarStr,'.','p'); % Replace period to avoid clash with Windows extensions

% Specify output directory and filename
outputDir = '.\Output\';
outputFilename = ['CMSgiven' compStr '_' TstarStr '.mat'];

% Load ERF and GMPM data for disaggregation
PSHAdir = '..\PSHA\Output\';
load(fullfile(PSHAdir, 'OpenSHAout.mat'));
load(fullfile(PSHAdir, 'GMPMout.mat'));

% Assumptions for ERF
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
Z2p5 = siteInfo.z2pt5;

% Load UHS data
load(fullfile(outputDir, 'uhs.mat'));
UHScurr_H = UHS_horiz_geo;
UHScurr_V = UHS_vert;

%% Define periods for CMS
Tcms = Tuhs; % Assume vibration periods specified for hazard curves include T*
nT = length(Tcms);
if ismember(Tstar,Tuhs)
    idTstar = find(Tuhs == Tstar);
else
    [~,id] = min(abs(Tuhs - Tstar));    
    TstarClosest = Tuhs(id);
    error(['Hazard curve period closest to T* is: ' num2str(TstarClosest) '; please revise Tstar variable to this value.']);
end

%% Extract PSHA data for conditioning component of GM and T*
% Get target acceleration
Atgt = exp(interp1( log(Tuhs), log(UHScurr_V), log(Tstar) )); % Depends on conditioning component of GM

% Determine GMPM output for all rupture scenarios
MUs_Tstar = MUs_V_all(:,idTstar); % Depends on conditioning component of GM
SIGs_Tstar = SIGs_V_all(:,idTstar); % Depends on conditioning component of GM

%% Disaggregation at conditioning component of GM and T*
% Determine percent contribution of each rupture scenario
rupRates = siteSeismicity.rupList(:,3); % Annual rate of rupture occurrence
numerator = lognpdf(Atgt,MUs_Tstar,SIGs_Tstar) .* rupRates; % Column vector
denominator = rupRates' * lognpdf(Atgt,MUs_Tstar,SIGs_Tstar); % Scalar
rupPercContr = numerator/denominator;

% Determine mean controlling scenario
meanScenario = rupPercContr' * siteSeismicity.rupList(:,[4 5]); % Mean magnitude and Rrup

% Display output
Mbar = meanScenario(1,1); Rbar = meanScenario(1,2);
fprintf('Mean scenario: M=%2.1f at R=%2.1f km\n',Mbar,Rbar);

%% Compute CMS - GMPM output for mean scenario
% GMPM output for H
AmH = zeros(1,nT); sH = zeros(1,nT); rhoHHTstarT = zeros(1,nT);
for ii=1:nT
    [AmH(1,ii), sH(1,ii)] = CB_2014_nga(Mbar, Tcms(ii), Rbar, Rbar, Rbar, W, Ztor, Zbot, dip, lambda, Fhw, Vs30, Z2p5, Zhyp, region);
    rhoHHTstarT(1,ii) = baker_jayaram_correlation(Tstar,Tcms(ii));    
end

% GMPM output for V/H
VHratio = zeros(1,nT); sVH = zeros(1,nT); 
rhoHVHTTstar = zeros(1,nT); rhoHVHsameT = zeros(1,nT);
for ii=1:nT
    [VHratio(1,ii), sVH(1,ii)] = gmpmVoverH_BC2016(Mbar,Rbar,Rbar,Rbar,lambda,dip,Fhw,Vs30,region,Sj,Tcms(ii),W,Ztor,Z2p5,Zhyp);
    rhoHVHTTstar(1,ii) = GA_2011_corr(Tcms(ii), Tstar, Mbar, Rbar, Vs30, FRV, FNM); % Note: rho_H_V/H(T,T*) not the same as rho_H_V/H(T*,T)
    rhoHVHsameT(1,ii) = GA_2011_corr(Tcms(ii), Tcms(ii), Mbar, Rbar, Vs30, FRV, FNM);
end

% GMPM output for V
AmV = zeros(1,nT); sV = zeros(1,nT); 
rhoVVTstarT = zeros(1,nT); 
for ii=1:nT
    [AmV(1,ii), sV(1,ii)] = gmpmV_BC2016(Mbar, Rbar, Rbar, Rbar, FRV, FNM, dip, Vs30, region, Sj, Tcms(ii), W, Ztor, Z2p5, Zhyp);    
    rhoVVTstarT(1,ii) = GKAS_2016_corr(Tstar, Tcms(ii), Mbar, Rbar, Rbar, Rbar, FRV, FNM, dip, Vs30, region, Sj); % Needed only for CMS given Av(T*)
end

% Derive data for V using H and V/H
AmV_derived = AmH.*VHratio;
sV_derived = sqrt( sH.^2 + sVH.^2 + 2*rhoHVHsameT.*sH.*sVH );
rhoHVTTstar_derived = (sH(idTstar).*rhoHHTstarT + sVH(idTstar).*rhoHVHTTstar) ./ sV_derived(idTstar); % Needed only for CMS given Av(T*)

%% Compute CMS - Epsilon at conditioning component of GM and T*
% Epsilon at T*
Am_Tstar = AmV(1,idTstar); % Depends on conditioning component of GM
sig_Tstar = sV(1,idTstar); % Depends on conditioning component of GM
epsTstar = (log(Atgt)-log(Am_Tstar)) / sig_Tstar;

%% Compute CMS - target median and dispersion
% CMS_H
CMS_H = AmH .* exp( sH .* (epsTstar * rhoHVTTstar_derived) ); % Depends on conditioning component of GM

% CMS_V
CMS_V = AmV .* exp( sV .* (epsTstar * rhoVVTstarT) ); % Depends on conditioning component of GM

% Conditional standard deviations
CMSsig_H = sH.*sqrt(1-rhoHVTTstar_derived.^2); % Depends on conditioning component of GM
CMSsig_V = sV.*sqrt(1-rhoVVTstarT.^2); % Depends on conditioning component of GM

%% Save data
save(fullfile(outputDir,outputFilename),...
    'Tstar','Tcms',...
    'Atgt','meanScenario',...
    'AmH','sH','rhoHHTstarT','VHratio','sVH','rhoHVHsameT','AmV','sV',...
    'rhoHVHTTstar','rhoVVTstarT','rhoHVTTstar_derived',...
    'epsTstar','CMS_H','CMS_V','CMSsig_H','CMSsig_V');