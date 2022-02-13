function [ IMs ] = find_ground_motions_wVert( selectionParams, targetSa, IMs )
%{
By: Neal Simon Kwong, nealsimonkwong@berkeley.edu

Paper citation:

Kwong, N.S. and Chopra, A.K. (2019). Selecting, scaling, and orienting
three components of ground motions for intensity-based assessments at
far-field sites. Earthquake Spectra. (in review).
%}

%% Additional inputs to select both H and V components of GM
% Define target spectrum for each GM in final selected ensemble
targetSpectra_H = repmat( exp( targetSa.meanReq      ), selectionParams.nGM, 1); % Single target spectrum for all GMs in final selected ensemble 
targetSpectra_V = repmat( exp( targetSa.meanReq_vert ), selectionParams.nGM, 1); % Single target spectrum for all GMs in final selected ensemble

% Relative importance of components of GM
wH = selectionParams.wH;
wV = 1-wH; 

%% Identify vibration periods for scaling GMs
if strcmp( selectionParams.condComp, 'H' )
    scaleFacIndex = selectionParams.indTcond; % Scale H spectrum to match, for example, T*
    scaleFacIndex_vert = (1:length(selectionParams.TgtPer))'; % Scale V spectrum to match full range of periods defined by selectionParams.TgtPer
else    
    scaleFacIndex_vert = selectionParams.indTcond; % Scale V spectrum to match, for example, T*
    scaleFacIndex = (1:length(selectionParams.TgtPer))'; % Scale H spectrum to match full range of periods defined by selectionParams.TgtPer
end

%% Select GMs
% Initialize
IMs.recID = zeros(selectionParams.nGM,1); % IDs of filtered database corresponding to selected GMs
IMs.sampleSmall = []; % Logarithmic spectra for H components of selected GMs at requested periods
IMs.sampleSmall_vert = []; % Logarithmic spectra for V components of selected GMs at requested periods
IMs.scaleFac = ones(selectionParams.nGM,1); % Scale factors for H components of selected GMs
IMs.scaleFac_vert = ones(selectionParams.nGM,1); % Scale factors for V components of selected GMs

% For each GM in final selected ensemble, find unique record in GM database that is closest to target spectrum
for i = 1:selectionParams.nGM
    err = zeros(selectionParams.nBig,1); % initialize error vector (SSD_combined)
    scaleFac = ones(selectionParams.nBig,1); % initialize scale factors to 1 (these are never changed if no scaling is allowed)
    scaleFac_vert = ones(selectionParams.nBig,1); % initialize scale factors to 1 (these are never changed if no scaling is allowed)

    % compute scale factors and errors for each candidate ground motion
    for j=1:selectionParams.nBig
        if selectionParams.isScaled % if scaling is allowed, compute scale factor (otherwise it is already defined as equal to 1)
            if strcmp( selectionParams.condComp, 'H' )
                scaleFac(j) = geomean( targetSpectra_H(i,scaleFacIndex) ./ exp(IMs.sampleBig(j,scaleFacIndex)) );
                if selectionParams.sameSF == 1
                    scaleFac_vert(j) = scaleFac(j);
                else
                    scaleFac_vert(j) = geomean( targetSpectra_V(i,scaleFacIndex_vert) ./ exp(IMs.sampleBig_vert(j,scaleFacIndex_vert)) ); % Can change this here to scale V component to different range of periods
                end
            else
                scaleFac_vert(j) = geomean( targetSpectra_V(i,scaleFacIndex_vert) ./ exp(IMs.sampleBig_vert(j,scaleFacIndex_vert)) );                
                if selectionParams.sameSF == 1
                    scaleFac(j) = scaleFac_vert(j);
                else
                    scaleFac(j) = geomean( targetSpectra_H(i,scaleFacIndex) ./ exp(IMs.sampleBig(j,scaleFacIndex)) ); % Can change this here to scale V component to different range of periods
                end
            end                      
        end        
        SSD_H = sum(( log(scaleFac(j)*(exp(IMs.sampleBig(j,:)))) - log(targetSpectra_H(i,:)) ).^2);
        SSD_V = sum(( log(scaleFac_vert(j)*(exp(IMs.sampleBig_vert(j,:)))) - log(targetSpectra_V(i,:) )).^2);
        err(j) = wH*SSD_H + wV*SSD_V;
    end
    
    % revise errors to deliberately exclude ground motions
    err(IMs.recID(1:i-1)) = 1e6; % exclude previously-selected ground motions using large error
    err(scaleFac>selectionParams.maxScale | scaleFac<(1/selectionParams.maxScale)) = 1e6; % exclude ground motions requiring too much scaling
    err(scaleFac_vert>selectionParams.maxScale | scaleFac_vert<(1/selectionParams.maxScale)) = 1e6; % exclude ground motions requiring too much scaling
        
    % find minimum-error ground motion    
    [tmp, IMs.recID(i)] = min(err);
    assert(tmp < 1e6, 'Warning: problem with simulated spectrum. No good matches found. If selectionParams.maxScale is set to 1, make sure selectionParams.isScaled is set to 0.');
    IMs.scaleFac(i) = scaleFac(IMs.recID(i)); % store scale factor
    IMs.sampleSmall = [IMs.sampleSmall; log(exp(IMs.sampleBig(IMs.recID(i),:))*scaleFac(IMs.recID(i)))]; % store scaled log spectrum
    IMs.scaleFac_vert(i) = scaleFac_vert(IMs.recID(i));
    IMs.sampleSmall_vert = [IMs.sampleSmall_vert; log(exp(IMs.sampleBig_vert(IMs.recID(i),:))*scaleFac_vert(IMs.recID(i)))]; % store scaled log spectrum    
end
end