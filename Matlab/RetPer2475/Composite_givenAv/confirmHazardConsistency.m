%% Confirm hazard consistency of selected GMs for all three components of GM
clear all; close all; clc;

%% Define case to be considered
GMSMinput

% Create labels
TminStr = num2str(Tmin); TminStr = strrep(TminStr,'.','p'); % Replace period to avoid clash with Windows extensions
TmaxStr = num2str(Tmax); TmaxStr = strrep(TmaxStr,'.','p'); % Replace period to avoid clash with Windows extensions
wHstr = num2str(wH); wHstr = strrep(wHstr,'.','p'); % Replace period to avoid clash with Windows extensions
caseStr = ['T' TminStr '_' TmaxStr '_SFmax' int2str(SFlimit) '_sameSF' int2str(sameSFflag) '_wH' wHstr];

%% Input
% Input directories and filenames
targetSpectraDir = '..\Output\';
targetSpectraMatFn = ['CompositeGiven' compStr '_' TminStr '_' TmaxStr '.mat'];
GMSMoutDir = ['.\GMSMout\' caseStr];
GMSMoutMatFn = 'GMSMout.mat';
timeSeriesOutMatFn = 'timeSeriesData.mat';
rotatedTimeSeriesMatFn = 'GMSMout_rotated.mat';

% Output directory
figDir = ['.\Figures\' caseStr];

%% Load target spectrum
load(fullfile(targetSpectraDir,targetSpectraMatFn),...
    'Tcms','CompSpec_H','CompSpec_V');
Ttgt = Tcms;
tgtSpec_H = CompSpec_H;
tgtSpec_V = CompSpec_V;

% Define bounds for comparing spectra
accLimit = 2; % Limit for response spectra
tgtBoundsLo_H = CompSpec_H ./ accLimit;
tgtBoundsHi_H = CompSpec_H .* accLimit;
tgtBoundsLo_V = CompSpec_V ./ accLimit;
tgtBoundsHi_V = CompSpec_V .* accLimit;

%% Load time series of selected GMs
% GMSM output (no time series)
load(fullfile(GMSMoutDir,GMSMoutMatFn),...
    'NGAselected','sfSelected','sfSelected_vert','metadata');

% GMSM output (yes time series, unrotated)
load(fullfile(GMSMoutDir,timeSeriesOutMatFn),...
    'gacc_all','dt_all','npts_all');
nGM = size(gacc_all,1); nComp = size(gacc_all,2);

% GMSM output (yes time series, rotated)
load(fullfile(GMSMoutDir,rotatedTimeSeriesMatFn),...
    'thetaOpt_all','gaccSubRot_all','npts_min');

% Replace time series for H components with those corresponding to optimal azimuths
gaccCurr = gacc_all; % The variable gaccCurr stores unmodified time series for all 3 components of GM for all selected records
gaccCurr(:,1:2) = gaccSubRot_all; 
    
%% Compute response spectra for scaled time series of selected GMs
% Initialize
T = metadata.knownPer; % Vibration periods for computing spectra
dampR = 0.05; % Damping ratio for computing spectra, unitless
nT = length(T);
scaledSpectra_all = zeros(nGM,nT,nComp); % Given component of GM, response spectra for all selected GMs

% Compute response spectra
for idComp = 1:nComp
    for idGM = 1:nGM      
        % Get unscaled time series
        accel_unscaled = gaccCurr{idGM,idComp};
        
        % Scale time series
        if idComp==3
            accel_scaled = accel_unscaled * sfSelected_vert(idGM,1);
        else
            accel_scaled = accel_unscaled * sfSelected(idGM,1);
        end
        
        % Get time step for time series
        dt_curr = dt_all(idGM,idComp);    
        
        % Compute response spectrum
        scaledSpectra_all(idGM,:,idComp) = calcResponseSpectrum(accel_scaled,dt_curr,T,dampR)';
    end
end       

%% Compare spectra
% Fig inputs
xMin = 1e-2; xMax = 1e1; yMin = 1e-3; yMax = 4e0;

% Plot
hFig = figure('units','inches','position',[0.1 0.5 10 6]);

% Horizontal 1
idComp = 1;
subplot(2,2,idComp);
specSelCurr = scaledSpectra_all(:,:,idComp);
medSpecSelCurr = geomean(specSelCurr,1);

hSel_H = loglog(T,specSelCurr,'color',ones(1,3)*.75);
hold on;
hTgt_H = loglog(Ttgt,tgtSpec_H,'color','k','linestyle','-','linewidth',2);
hAccLimit = plot(Ttgt,tgtBoundsLo_H,'color','k','linestyle','--','linewidth',2);
plot(Ttgt,tgtBoundsHi_H,'color','k','linestyle','--','linewidth',2);
hSelMed_H = plot(T,medSpecSelCurr,'r-.','linewidth',2);

% Guideline
line(ones(1,2)*Tmin,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmin,yMax,'\itT\rm_{min}','verticalalignment','top','horizontalalignment','right');
line(ones(1,2)*Tmax,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmax,yMax,'\itT\rm_{max}','verticalalignment','top','horizontalalignment','left');

% Format axes
axis([xMin xMax yMin yMax]);

% Labels
ylabel('Spectral acceleration, g'); 
title('H1 component');
hLeg = legend([hTgt_H hAccLimit hSel_H(1) hSelMed_H],...
    {'Composite Spectrum' 'Acceptable limits' 'Selected GMs' 'Geomean of selected'});
set(hLeg,'box','off','location','southwest');


% Horizontal 2
idComp = 2;
subplot(2,2,idComp);
specSelCurr = scaledSpectra_all(:,:,idComp);
medSpecSelCurr = geomean(specSelCurr,1);

hSel_H = loglog(T,specSelCurr,'color',ones(1,3)*.75);
hold on;
hTgt_H = loglog(Ttgt,tgtSpec_H,'color','k','linestyle','-','linewidth',2);
hAccLimit = plot(Ttgt,tgtBoundsLo_H,'color','k','linestyle','--','linewidth',2);
plot(Ttgt,tgtBoundsHi_H,'color','k','linestyle','--','linewidth',2);
hSelMed_H = plot(T,medSpecSelCurr,'r-.','linewidth',2);

% Guideline
line(ones(1,2)*Tmin,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmin,yMax,'\itT\rm_{min}','verticalalignment','top','horizontalalignment','right');
line(ones(1,2)*Tmax,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmax,yMax,'\itT\rm_{max}','verticalalignment','top','horizontalalignment','left');

% Format axes
axis([xMin xMax yMin yMax]);

% Labels
xlabel('Period, sec'); title('H2 component');


% Vertical
idComp = 3;
subplot(2,2,idComp);
specSelCurr = scaledSpectra_all(:,:,idComp);
medSpecSelCurr = geomean(specSelCurr,1);

hSel_V = loglog(T,specSelCurr,'color',ones(1,3)*.75);
hold on;
hTgt_V = loglog(Ttgt,tgtSpec_V,'color','k','linestyle','-','linewidth',2);
hAccLimit = plot(Ttgt,tgtBoundsLo_V,'color','k','linestyle','--','linewidth',2);
plot(Ttgt,tgtBoundsHi_V,'color','k','linestyle','--','linewidth',2);
hSelMed_V = plot(T,medSpecSelCurr,'r-.','linewidth',2);

% Guideline
line(ones(1,2)*Tmin,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmin,yMax,'\itT\rm_{min}','verticalalignment','top','horizontalalignment','right');
line(ones(1,2)*Tmax,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmax,yMax,'\itT\rm_{max}','verticalalignment','top','horizontalalignment','left');

% Format axes
axis([xMin xMax yMin yMax]);

% Labels
xlabel('Period, sec'); ylabel('Spectral acceleration, g'); title('V component');

% Save fig
figFn = ['hazConsConfirmation']; 
figFnWithPath_pdf = fullfile(figDir,[figFn '.pdf']);
set(hFig,'color','w'); export_fig(figFnWithPath_pdf);

%% Summarize GMSM
GMSMtable = [NGAselected sfSelected sfSelected_vert thetaOpt_all];
open GMSMtable