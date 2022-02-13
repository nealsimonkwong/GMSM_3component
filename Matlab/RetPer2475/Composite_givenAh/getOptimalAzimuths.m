%% Determine optimal azimuths for H components of selected records
clear all; close all; clc;

%% Define case to be considered
GMSMinput

% Create labels
TminStr = num2str(Tmin); TminStr = strrep(TminStr,'.','p'); % Replace period to avoid clash with Windows extensions
TmaxStr = num2str(Tmax); TmaxStr = strrep(TmaxStr,'.','p'); % Replace period to avoid clash with Windows extensions
wHstr = num2str(wH); wHstr = strrep(wHstr,'.','p'); % Replace period to avoid clash with Windows extensions
caseStr = ['T' TminStr '_' TmaxStr '_SFmax' int2str(SFlimit) '_sameSF' int2str(sameSFflag) '_wH' wHstr];

%% Inputs
% Input directories and filenames
targetSpectraDir = '..\Output\';
GMSMoutDir = ['.\GMSMout\' caseStr];
GMSMoutMatFn = 'GMSMout.mat';
timeSeriesOutMatFn = 'timeSeriesData.mat';

% Output filename
rotatedTimeSeriesMatFn = 'GMSMout_rotated.mat';

% Load target spectrum
outputFilename = ['CompositeGiven' compStr '_' TminStr '_' TmaxStr '.mat'];
load(fullfile(targetSpectraDir,outputFilename),...
    'Tcms','CompSpec_H','CompSpec_V');
Ttgt = Tcms;
tgtSpectrum_x = CompSpec_H;
tgtSpectrum_y = CompSpec_H;

% Specify period range for comparing spectra
load(fullfile(GMSMoutDir,GMSMoutMatFn),'selectionParams');
nGM = selectionParams.nGM;
Tselection = selectionParams.TgtPer; % Same range and number of periods as those used when comparing spectra to select records

% Load GMSM output (as-recorded)
load(fullfile(GMSMoutDir,timeSeriesOutMatFn),...
    'gacc_all','dt_all','npts_all');
load(fullfile(GMSMoutDir,GMSMoutMatFn),'sfSelected');

% Inputs for computing response spectra
dampR = 0.05; % Damping ratio, unitless
thetaList = 0:1:179; % Specify range of azimuths in 1 deg increments
nTheta = length(thetaList);
nPer = length(Tselection);
nComp = 2;

%% Ensure number of data points in time series is identical for both H comp
% Check if dt is same for both H components
if abs(dt_all(:,1)-dt_all(:,2))>1e-12
    error('Time steps are different for each H component.\n');
end

% Use min npts among both H comps to truncate time series
npts_min = min(npts_all,[],2);
gaccSub_all = cell(nGM,nComp);
for ii=1:nGM
    for jj=1:nComp
        gaccSub_all{ii,jj} = gacc_all{ii,jj}(1:npts_min(ii,1));
    end
end

%% Compute response spectra for all orientations of each GM pair
% Given GM and vibration period, compute pseudo-acceleration histories of
% SDF system due to different orientations of input GM using superposition
tic
% Initialize
rotatedSpectra_all = cell(nGM,nComp); % Given GM record and H component, store response spectra corresponding to all non-redundant angles (nTheta x nPer matrix)
for idGM = 1:nGM
    % Fix scaled GM pair
    uddgx = gaccSub_all{idGM,1} * sfSelected(idGM); % units of g
    uddgy = gaccSub_all{idGM,2} * sfSelected(idGM); % units of g
    dt = dt_all(idGM,1); % dt should be same for both H components
        
    % Initialize storage vars for spectra corresponding to all rotated versions of GM pair
    currRotSpec_x = zeros(nTheta,nPer);
    currRotSpec_y = zeros(nTheta,nPer);
    % Given vibration period, apply superposition to determine
    % pseudo-acceleration histories for all non-redundant angles
    for ii = 1:nPer
        % Find as-recorded pseudo-acceleration histories once
        Tcurr = Tselection(ii);
        wn = 2*pi/Tcurr;
        AhistX = LinearInterpolation(1,wn,dampR,0,0,uddgx,dt,1)' * wn^2;
        AhistY = LinearInterpolation(1,wn,dampR,0,0,uddgy,dt,1)' * wn^2;        
        % Derive pseudo-acceleration histories for rotated GM pairs via superposition
        for tt = 1:nTheta
            theta = thetaList(tt);
            AhistRot = AhistX*cosd(theta) - AhistY*sind(theta); % Superposition
            currRotSpec_x(tt,ii) = max(abs( AhistRot )); % Get spectral ordinate
        end
    end
    
    % Use periodicity to complete spectra for orthogonal H component
    currRotSpec_y( thetaList<90 ,:) = currRotSpec_x( thetaList>=90 ,:);
    currRotSpec_y( thetaList>=90 ,:) = currRotSpec_x( thetaList<90 ,:);
    
    % Save for given GM pair
    rotatedSpectra_all{idGM,1} = currRotSpec_x;
    rotatedSpectra_all{idGM,2} = currRotSpec_y;
end
toc

%% Save spectra for all rotated versions of GM pairs
save(fullfile(GMSMoutDir,rotatedTimeSeriesMatFn),...
    'dampR','Tselection','thetaList',...
    'rotatedSpectra_all'); 

%% Determine optimal azimuth for each GM using pair of target spectra
% Put pair of target spectra under common period scale
tgtSpectrum_x_sub = exp( interp1(log(Ttgt),log(tgtSpectrum_x),log(Tselection)) );
tgtSpectrum_y_sub = exp( interp1(log(Ttgt),log(tgtSpectrum_y),log(Tselection)) );

% Initialize vars
thetaOpt_all = zeros(nGM,1); % List of angles for rotating (CCW from plan view) each GM pair
spectraSelectedRot_x = zeros(nGM,nPer);
spectraSelectedRot_y = zeros(nGM,nPer);

% Fix GM pair
for idGM = 1:nGM
    % Get all rotated spectra for current GM pair
    spectra_x_sub = rotatedSpectra_all{idGM,1};
    spectra_y_sub = rotatedSpectra_all{idGM,2};
    
    % Compute misfits for all rotated spectra
    SSD_x = sum((   bsxfun(@minus, log(tgtSpectrum_x_sub), log(spectra_x_sub)  )   ).^2,2);
    SSD_y = sum((   bsxfun(@minus, log(tgtSpectrum_y_sub), log(spectra_y_sub)  )   ).^2,2);
    SSD = SSD_x + SSD_y;
    
    % Find best azimuth with minimal misfit
    [~,sortedIndices] = sort(SSD);
    thetaOpt_all(idGM,1) = sortedIndices(1);
    spectraSelectedRot_x(idGM,:) = spectra_x_sub(sortedIndices(1),:);
    spectraSelectedRot_y(idGM,:) = spectra_y_sub(sortedIndices(1),:);
end

%% Save optimal azimuths
save(fullfile(GMSMoutDir,rotatedTimeSeriesMatFn),...
    'thetaOpt_all',...
    'spectraSelectedRot_x','spectraSelectedRot_y','-append');

%% Determine unscaled GM pairs rotated CCW by angles corresponding to optimal azimuths
gaccSubRot_all = cell(nGM,nComp);
for idGM = 1:nGM
    % Fix scaled GM pair
    a1 = gaccSub_all{idGM,1}; % units of g
    a2 = gaccSub_all{idGM,2}; % units of g
    
    % Find rotated GM pair
    theta = thetaOpt_all(idGM,1);
    aRot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)] * [a1 a2]'; % Rotate scaled GM pair CCW by theta
    gaccSubRot_all{idGM,1} = aRot(1,:)';
    gaccSubRot_all{idGM,2} = aRot(2,:)';     
end

%% Save rotated GM pairs
save(fullfile(GMSMoutDir,rotatedTimeSeriesMatFn),...
    'npts_min',...
    'gaccSub_all','gaccSubRot_all','-append');