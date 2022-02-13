function [ selectionParams, metadata, indPer ] = screen_database(selectionParams, allowedRecs )
% Load a database of ground motion data, and screen it to identify usable
% ground motions for potential selection

% load the specified database
load(fullfile(selectionParams.databaseDir,selectionParams.databaseFile));

% Format appropriate ground motion metadata variables. Additional metadata 
% is available in the databases and can be added here if desired.
% Note: These lines should be modified if using the BBP_EXSIM_meta_data.mat
% database file. See documentation for more details.
metadata.getTimeSeries = getTimeSeries; % Instructions for getting time series from NGAW2
metadata.recNum      = [1:length(NGA_num)]; % Total number of records in database
if selectionParams.RotD == 50 && exist('Sa_RotD50')
    SaKnown     = Sa_RotD50;
elseif selectionParams.RotD == 100 && exist('Sa_RotD100')
    SaKnown     = Sa_RotD100;
else
    display(['Error--RotD' num2str(RotD) ' not provided in database'])
    % If data corresponding to user input RotD value does not exist,
    % optionally use the geometric mean of the single-component Sa's:
    % SaKnown = sqrt(Sa_1.*Sa_2);
end
SaKnown_vert = Sa_vert;

%% Arrange available spectra in usable format and check for invalid values
% Create variable for known periods
TmaxGMPM = 10; % 10 or 4 sec depending on GMPM
idxPer = find(Periods <= TmaxGMPM); % Exclude periods that cause problems with evaluating a GMPM
knownPer = Periods(idxPer); 

% Ensure specified T* is an element of knownPer
if ~ismember(selectionParams.Tcond,knownPer)
    disp('Warning: Conditioning period is not one of the periods used for pre-computing GM spectra');
end

% Modify TgtPer to include Tcond if running a conditional selection
if selectionParams.cond == 1 && any(ismember(selectionParams.Tcond,selectionParams.TgtPer)==0) % If any Tstar not in TgtPer  
    selectionParams.TgtPer = sort([selectionParams.TgtPer selectionParams.Tcond]); % Include
end

% Match periods (known periods and target periods for SSD computations) 
% save the indices of the matched periods in knownPer
indPer = zeros(length(selectionParams.TgtPer),1);
for i=1:length(selectionParams.TgtPer)
    [~ , indPer(i)] = min(abs(knownPer - selectionParams.TgtPer(i)));
end

% Remove any repeated values from TgtPer and update TgtPer using periods 
% provided in GM database
indPer = unique(indPer);
selectionParams.TgtPer = knownPer(indPer);

% Identify the indices of Tcond within TgtPer (used for scaling GMs later)
temp = zeros(size(selectionParams.Tcond));
for ii=1:length(selectionParams.Tcond)
    [~, temp(ii)] = min(abs(selectionParams.TgtPer - selectionParams.Tcond(ii)));
end
selectionParams.indTcond = temp;

%% Screen the records to be considered
% Determine valid records within GM database
recValidSa = ~any(Sa_1==-999,2) & ~any(Sa_2==-999,2) & ~any(ismember(Sa_vert,[-999 inf]),2); % Ensure that each of remaining records contains all 3 components
recValidSoil = soil_Vs30 > allowedRecs.Vs30(1) & soil_Vs30 < allowedRecs.Vs30(2);
recValidMag =  magnitude > allowedRecs.Mag(1)  & magnitude < allowedRecs.Mag(2);
recValidDist = closest_D > allowedRecs.D(1)    & closest_D < allowedRecs.D(2);
recValidLUF = lowest_usable_freq > allowedRecs.LUF(1) & lowest_usable_freq < allowedRecs.LUF(2);
recValidNGA = ~ismember(NGA_num,allowedRecs.NGAinvalid);
metadata.allowedIndex = find(recValidSa & recValidSoil & recValidMag & recValidDist & recValidLUF & recValidNGA); 

% resize SaKnown to include only allowed records
SaKnown = SaKnown(metadata.allowedIndex,idxPer);       
SaKnown_vert = SaKnown_vert(metadata.allowedIndex,idxPer); % Assume spectra for V component computed at same periods as spectra for H component

% count number of allowed spectra in reduced GM database
selectionParams.nBig = length(metadata.allowedIndex);  
display(['Number of remaining ground motions for selection = ' num2str(selectionParams.nBig)])
assert(selectionParams.nBig >= selectionParams.nGM, 'Warning: there are not enough allowable ground motions');

% Save spectra of remaining GMs in the database
metadata.knownPer = knownPer;
metadata.SaKnown = SaKnown;
metadata.SaKnown_vert = SaKnown_vert;
end