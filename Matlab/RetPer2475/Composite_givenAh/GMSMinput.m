%% Specify inputs for using scripts under this "target spectrum" folder
% Target spectrum
tgtSpecStr = 'Composite';
compStr = 'Ah'; % Which component to condition upon
compForConditioning = 'H'; % Which component to condition upon

% Conditioning periods
Tmin = 0.1; 
Tmax = 0.75; 

% GMSM input
isScaledFlag = 1; % =1 scales GMs; otherwise, =0 and set SFlimit to unity
SFlimit = 4; % Limit scale factors (SFs) to within (1/SFlimit) to SFlimit
sameSFflag = 1; % =1 sets SF for V component to be same as SF for H component; otherwise, compute separate SFs for V component
wH = 0.5; % Ranges from 0 (only V component) to 1 (only H component) when comparing spectra using SSDcombined for GM selection