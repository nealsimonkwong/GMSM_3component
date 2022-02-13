%{ 
Function to select ground motions from a given set of inputs
Outputs variables that are updated by function (selectionParams and
targetSa)
%}
function [IMs, selectionParams, targetSa, metadata] = selectGMs_MAIN( selectionParams, allowedRecs, targetSa)
% Load the specified ground motion database and screen for suitable motions
[selectionParams, metadata, indPer] = screen_database( selectionParams, allowedRecs );
IMs.sampleBig = log(metadata.SaKnown(:,indPer));  % store logarithmic spectral accelerations at vibration periods closest to requested periods (selectionParams.TgtPer)
IMs.sampleBig_vert = log(metadata.SaKnown_vert(:,indPer));

% Determine target spectrum at same periods as GM spectra
targetSa.meanReq = interp1( log(targetSa.Ttgt), log(targetSa.tgtSpec_H), log(metadata.knownPer(indPer)) );
targetSa.meanReq_vert = interp1( log(targetSa.Ttgt), log(targetSa.tgtSpec_V), log(metadata.knownPer(indPer)) );

% Find best matches to the simulated spectra from ground-motion database
IMs = find_ground_motions_wVert( selectionParams, targetSa, IMs );
end