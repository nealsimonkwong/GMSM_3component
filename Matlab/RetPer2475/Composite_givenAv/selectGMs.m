%% Select GMs to match a multicomponent target spectrum
% Matlab files from Prof. Jack Baker modified to include V component of GM
clear all; close all; clc;

%% Specify case to be considered
GMSMinput

% Create labels
TminStr = num2str(Tmin); TminStr = strrep(TminStr,'.','p'); % Replace period to avoid clash with Windows extensions
TmaxStr = num2str(Tmax); TmaxStr = strrep(TmaxStr,'.','p'); % Replace period to avoid clash with Windows extensions
wHstr = num2str(wH); wHstr = strrep(wHstr,'.','p'); % Replace period to avoid clash with Windows extensions
caseStr = ['T' TminStr '_' TmaxStr '_SFmax' int2str(SFlimit) '_sameSF' int2str(sameSFflag) '_wH' wHstr]; % Create unique names for different cases of GM selection

% Input directories
targetSpectraDir = '..\Output\';

% Output directories
outDir = ['.\GMSMout\' caseStr];
figDir = ['.\Figures\' caseStr];
if(~exist(outDir,'dir')); mkdir(outDir); end;
if(~exist(figDir,'dir')); mkdir(figDir); end;

% Output filename
outFn = ['GMSMout.mat'];
NGAidFn = ['NGAid.txt']; % Text file containing NGA Record Sequence Numbers (RSNs) of selected GMs


%% Load target spectrum
outputFilename = ['CompositeGiven' compStr '_' TminStr '_' TmaxStr '.mat'];
load(fullfile(targetSpectraDir,outputFilename),...
    'Tcms','CompSpec_H','CompSpec_V');
targetSa.Ttgt = Tcms;
targetSa.tgtSpec_H = CompSpec_H;
targetSa.tgtSpec_V = CompSpec_V;

% Define bounds for comparing spectra
accLimit = 2; % Limit for response spectra (with sigma approximately 0.7, exponent of sigma is approximately two)
targetSa.tgtBoundsLo_H = CompSpec_H ./ accLimit;
targetSa.tgtBoundsHi_H = CompSpec_H .* accLimit;
targetSa.tgtBoundsLo_V = CompSpec_V ./ accLimit;
targetSa.tgtBoundsHi_V = CompSpec_V .* accLimit;


%% Specify GMSM inputs
% Ground motion database and type of selection 
selectionParams.databaseDir     = '..\..\Databases';
selectionParams.databaseFile    = 'NGA_W2_meta_data';
selectionParams.RotD            = 50; % Define H component of GM
selectionParams.nGM             = 11; % number of ground motions to be selected 

% Spectral periods of interest
selectionParams.cond       = 1; % =1 to include T* into periods for selecting GMs (selectionParams.TgtPer)
selectionParams.Tcond      = logspace(log10(Tmin),log10(Tmax),20); % Periods where spectra should be "matched" via scaling; e.g., T* or perhaps a wider range like Tmin to Tmax
selectionParams.TgtPer     = logspace(log10(0.01),log10(10),20); % Periods for comparing spectra in selection; this variable will be updated with periods where GM spectra already computed that are also closest to requested periods; see screen_database() function

% Parameters for scaling GMs
selectionParams.condComp   = compForConditioning; % ='H' to scale GMs to match H component or ='V' to scale GMs to match V component
selectionParams.isScaled   = isScaledFlag; % =1 scales GMs; otherwise, =0 and set SFlimit to unity
selectionParams.maxScale   = SFlimit; % Limit scale factors (SFs) to within (1/SFlimit) to SFlimit
selectionParams.sameSF     = sameSFflag; % =1 sets SF for V component to be same as SF for H component; otherwise, compute separate SFs for V component

% Relative importance of H component of GM
selectionParams.wH         = wH; % Ranges from 0 (only V component) to 1 (only H component) when comparing spectra using SSDcombined for GM selection

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
% Limits on metadata for filtering GM database
allowedRecs.Mag  = [ 5  Inf];       % upper and lower bound of allowable magnitude values
allowedRecs.Vs30 = [180 Inf];       % upper and lower bound of allowable Vs30 values 
allowedRecs.LUF  = [ 0  0.1];       % upper and lower bound of allowable Lowest usable frequency values
allowedRecs.D    = [ 15  Inf];      % upper and lower bound of allowable distance values
allowedRecs.NGAinvalid = [4577:4839 6993:8055 9194]; % Exclude NGA Record Sequence Numbers (RSNs) from consideration for various reasons (e.g., bc cannot retrieve time series for these from PEER website)


%% Select GMs using functions
[IMs, selectionParams, targetSa, metadata] = selectGMs_MAIN( selectionParams, allowedRecs, targetSa);


%% Extract data from selection
% NGA RSNs
NGAselected = metadata.recNum( metadata.allowedIndex( IMs.recID ) )'; % RSNs of selected GMs, sorted from best to worst match relative to target spectrum

% Scale factors
sfSelected = IMs.scaleFac; % Scale factors for H components of selected GMs
sfSelected_vert = IMs.scaleFac_vert; % Scale factors for V components of selected GMs

% Determine spectra of selected (and scaled) GMs at periods given in GM
% database
spectraSelected = diag(IMs.scaleFac) * metadata.SaKnown(IMs.recID,:); % H components
spectraSelected_vert = diag(IMs.scaleFac_vert) * metadata.SaKnown_vert(IMs.recID,:); % V components

% Names of corresponding earthquakes
load(fullfile(selectionParams.databaseDir,selectionParams.databaseFile),'NGA_num','EQ_name');
idRec = find( ismember(NGA_num,NGAselected) );
EQselected = EQ_name(idRec);


%% Sort results based on ascending NGA IDs and save
% Sort
[~,ix] = sort(NGAselected,'ascend');
NGAselected = NGAselected(ix,1);
sfSelected = sfSelected(ix,1);
sfSelected_vert = sfSelected_vert(ix,1);
spectraSelected = spectraSelected(ix,:);
spectraSelected_vert = spectraSelected_vert(ix,:);
EQselected = EQselected(ix,:);

% Save mat file
save(fullfile(outDir,outFn),...
    'selectionParams','targetSa',...
    'IMs','metadata',...
    'NGAselected','sfSelected','sfSelected_vert','spectraSelected','spectraSelected_vert','EQselected');  

% Save selected NGA RSNs as text file for retrieving time series from PEER website
fid = fopen( fullfile(outDir,NGAidFn) ,'w');
for ii = 1:selectionParams.nGM
    fprintf(fid,'%i ,\n',NGAselected(ii));
end
fclose(fid);


%% Preliminary check before retrieving time series: Plot spectra of selected GMs against targets
% Fig inputs
xMin = 1e-2; xMax = 1e1; yMin = 1e-3; yMax = 2e0;

% Plot
hFig = figure('units','inches','position',[0.1 0.5 9 4.5]);
% Horizontal
subplot(1,2,1);
hSel_H = loglog(metadata.knownPer,spectraSelected,'color',ones(1,3)*.75);
hold on;
hTgtSpec = loglog(targetSa.Ttgt,targetSa.tgtSpec_H,'color','k','linestyle','-','linewidth',2);
hAccLimit = plot(targetSa.Ttgt,targetSa.tgtBoundsLo_H,'color','k','linestyle','--','linewidth',2);
plot(targetSa.Ttgt,targetSa.tgtBoundsHi_H,'color','k','linestyle','--','linewidth',2);
hSelMed_H = plot(metadata.knownPer,geomean(spectraSelected,1),'r-.','linewidth',2);

% Guideline
line(ones(1,2)*Tmin,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmin,yMin,'\itT\rm_{min}','verticalalignment','bottom','horizontalalignment','right');
line(ones(1,2)*Tmax,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmax,yMin,'\itT\rm_{max}','verticalalignment','bottom','horizontalalignment','left');

% Format axes
axis square;
axis([xMin xMax yMin yMax]);

% Labels
xlabel('Period, sec'); ylabel('Spectral acceleration, g'); title('Horizontal');
hLeg = legend([hTgtSpec hAccLimit hSel_H(1) hSelMed_H],...
    {'Target spectrum' 'Acceptable limits' 'Selected GMs' 'Geomean of selected'});
set(hLeg,'box','off','location','west');


% Vertical
subplot(1,2,2);
hSel_V = loglog(metadata.knownPer,spectraSelected_vert,'color',ones(1,3)*.75);
hold on;
hCompSpec_V = loglog(targetSa.Ttgt,targetSa.tgtSpec_V,'color','k','linestyle','-','linewidth',2);
hAccLimit = plot(targetSa.Ttgt,targetSa.tgtBoundsLo_V,'color','k','linestyle','--','linewidth',2);
plot(targetSa.Ttgt,targetSa.tgtBoundsHi_V,'color','k','linestyle','--','linewidth',2);
hSelMed_V = plot(metadata.knownPer,geomean(spectraSelected_vert,1),'r-.','linewidth',2);

% Guideline
line(ones(1,2)*Tmin,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmin,yMin,'\itT\rm_{min}','verticalalignment','bottom','horizontalalignment','right');
line(ones(1,2)*Tmax,[yMin yMax],'color','k','linestyle','-.'); 
text(Tmax,yMin,'\itT\rm_{max}','verticalalignment','bottom','horizontalalignment','left');

% Format axes
axis square;
axis([xMin xMax yMin yMax]);

% Labels
xlabel('Period, sec'); title('Vertical');

% Save fig
figFn = ['hazConsCheck']; 
figFnWithPath_pdf = fullfile(figDir,[figFn '.pdf']);
set(hFig,'color','w'); export_fig(figFnWithPath_pdf);