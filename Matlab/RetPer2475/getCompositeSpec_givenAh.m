%% Get multicomponent Composite Spectrum given Ah at Tmin to Tmax
% Note: Ensure beforehand that CMSs have been determined for Tmin and Tmax
clear all; close all; clc;

%% Input
% Define component of GM for conditioning
compStr = 'Ah';

% Specify range of vibration periods
Tmin = 0.1; % Make sure CMS has been computed for this period
Tmax = 0.75; % Make sure CMS has been computed for this period

% Create label for chosen Tmin and Tmax
TminStr = num2str(Tmin);
TminStr = strrep(TminStr,'.','p'); % Replace period to avoid clash with Windows extensions
TmaxStr = num2str(Tmax);
TmaxStr = strrep(TmaxStr,'.','p'); % Replace period to avoid clash with Windows extensions

% Specify output directory and filename
outputDir = '.\Output\';
outputFilename = ['CompositeGiven' compStr '_' TminStr '_' TmaxStr '.mat'];

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

%% Load CMSs for Tmin and Tmax
% Initialize
nT = length(Tuhs);
CMSs_H = zeros(2,nT);
CMSs_V = zeros(2,nT);

% Load CMS output for Tmin
load(fullfile(outputDir, ['CMSgiven' compStr '_' TminStr '.mat']));
CMSs_H(1,:) = CMS_H;
CMSs_V(1,:) = CMS_V;

% Load CMS output for Tmax
load(fullfile(outputDir, ['CMSgiven' compStr '_' TmaxStr '.mat']));
CMSs_H(2,:) = CMS_H;
CMSs_V(2,:) = CMS_V;

%% Conduct disaggregation at Tbar
% Determine Tbar
TbarExact = sqrt(Tmin*Tmax);
[~,idTbar] = min(abs(Tuhs - TbarExact));
Tbar = Tuhs(idTbar);
display(['Hazard curve period closest to Tbar=' num2str(TbarExact) ' is: ' num2str(Tbar)]);

% Get target acceleration
Atgt = UHScurr_H(1,idTbar); % Depends on conditioning component of GM

% Determine GMPM output for all rupture scenarios
MUs_Tbar = MUs_geo_all(:,idTbar); % Depends on conditioning component of GM
SIGs_Tbar = SIGs_geo_all(:,idTbar); % Depends on conditioning component of GM

% Determine percent contribution of each rupture scenario
rupRates = siteSeismicity.rupList(:,3); % Annual rate of rupture occurrence
numerator = lognpdf(Atgt,MUs_Tbar,SIGs_Tbar) .* rupRates; % Column vector
denominator = rupRates' * lognpdf(Atgt,MUs_Tbar,SIGs_Tbar); % Scalar
rupPercContr = numerator/denominator;

% Determine mean controlling scenario
meanScenario = rupPercContr' * siteSeismicity.rupList(:,[4 5]); % Mean magnitude and Rrup
Mbar = meanScenario(1,1); Rbar = meanScenario(1,2);

%% Get GMPM outputs for controlling mean scenario
% GMPM output for H        
[Ah, sigH] = CB_2014_nga(Mbar, Tcms, Rbar, Rbar, Rbar, W, Ztor, Zbot, dip, lambda, Fhw, Vs30, Z2p5, Zhyp, region);

% GMPM output for V
[Av, sigV] = gmpmV_BC2016(Mbar, Rbar, Rbar, Rbar, FRV, FNM, dip, Vs30, region, Sj, Tcms, W, Ztor, Z2p5, Zhyp);

% V/H output
VHratio = zeros(1,nT); sigVH = zeros(1,nT); 
rho_H_VH_sameT = zeros(1,nT);
for ii=1:nT
%     [VHratio(1,ii), sigVH(1,ii)] = GA_2011(Mbar, Rbar, Vs30, FRV, FNM, Tcms(1,ii));
    [VHratio(1,ii), sigVH(1,ii)] = gmpmVoverH_BC2016(Mbar,Rbar,Rbar,Rbar,lambda,dip,Fhw,Vs30,region,Sj,Tcms(1,ii),W,Ztor,Z2p5,Zhyp);            
    rho_H_VH_sameT(1,ii) = GA_2011_corr(Tcms(1,ii), Tcms(1,ii), Mbar, Rbar, Vs30, FRV, FNM);
end

% Derive rho_H_V_sameT
sigV_derived = sqrt( sigH.^2 + sigVH.^2 + 2*rho_H_VH_sameT.*sigH.*sigVH );
rho_H_V_sameT_derived = (sigH + sigVH.*rho_H_VH_sameT) ./ sigV_derived;

%% UHSbar (approximation using single disaggregation and single epsilon at Tbar)
% Note: For more accurate UHSbar, compute epsilons at every vibration
% period instead of only at Tbar
% Epsilon at Tbar
epsH = (log(UHScurr_H(idTbar) ./ Ah(idTbar))) ./ sigH(idTbar); % Depends on conditioning component of GM

% Target spectrum
UHSbarV = Av .* exp( epsH .* rho_H_V_sameT_derived .* sigV ); % Depends on conditioning component of GM

%% Composite Spectrum 
% Identify period ranges
idShort = find( Tcms<=Tmin );
idMed = find( Tcms>Tmin & Tcms<Tmax );
idLong = find( Tcms>=Tmax );

% H (RotD50)
CompSpec_H = zeros(1,nT);
CompSpec_H(idShort) = CMSs_H( 1 , idShort );
CompSpec_H(idMed) = UHScurr_H( 1 , idMed ); % Depends on conditioning component of GM
CompSpec_H(idLong) = CMSs_H( 2 , idLong );

% V (approx)
CompSpec_V = zeros(1,nT);
CompSpec_V(idShort) = CMSs_V( 1 , idShort );
CompSpec_V(idMed) = UHSbarV( 1 , idMed ); % Depends on conditioning component of GM
CompSpec_V(idLong) = CMSs_V( 2 , idLong );

%% Save data
save(fullfile(outputDir,outputFilename),...
    'Tmin','Tmax','Tcms',...
    'TbarExact','Tbar','UHSbarV',...
    'CompSpec_H','CompSpec_V');