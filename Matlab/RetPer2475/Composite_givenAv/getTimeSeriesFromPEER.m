%% Retrieve data from PEER time series
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
GMSMoutDir = ['.\GMSMout\' caseStr];
GMSMoutMatFn = 'GMSMout.mat';
PEERoutFn = 'PEERNGARecords_Unscaled';
NGAwest2outFn = '_SearchResults.csv';
RecordSetPath_orig = fullfile(GMSMoutDir,PEERoutFn); % Specify location of time series from PEER NGAWest2

% Output filename
timeSeriesOutMatFn = 'timeSeriesData.mat';

% Load GMSM output
load(fullfile(GMSMoutDir,GMSMoutMatFn));
NGAid_check = sort(NGAselected);

%% Gather info about all selected records from csv file
NGAid_csv = xlsread(fullfile(RecordSetPath_orig,NGAwest2outFn),'C35:C45'); % NGA RSNs
[~,at2fn_h1] = xlsread(fullfile(RecordSetPath_orig,NGAwest2outFn),'T35:T45'); % H1 filenames
[~,at2fn_h2] = xlsread(fullfile(RecordSetPath_orig,NGAwest2outFn),'U35:U45'); % H2 filenames
[~,at2fn_v] = xlsread(fullfile(RecordSetPath_orig,NGAwest2outFn),'V35:V45'); % V filenames
Ngm = length(NGAid_csv);

% Check if correct time series were obtained from NGAwest2
if all( NGAid_check == NGAid_csv ) ~= 1
    error('Wrong set of records obtained from PEER website.\n');
end

%% Read AT2 files for time series and extract GM info
% Horizontal 1
npts_all_h1 = zeros(Ngm,1); % Number of points in time series
dt_all_h1 = zeros(Ngm,1); % Time step between points in time series
gacc_all_h1 = cell(Ngm,1); % Ground acceleration data
for idNGA = 1:Ngm
    clear gacc
    currFn = at2fn_h1{idNGA}(2:end); % Str from csv file contains space in first entry
    fid = fopen(fullfile(RecordSetPath_orig,currFn),'r');
    for ll=1:3; fgets(fid); end; % Skip first 3 lines
    currLine = fgets(fid); % Get line containing data
    C = strsplit(currLine,', '); % Split line into two parts
    nptStr = C{1}; dtStr = C{2};
    temp = strsplit(nptStr,' '); % Further split string into parts
    npts = str2double( temp{2} );
    npts_all_h1(idNGA,1) = npts;
    temp = strsplit(dtStr,' '); % Further split string into parts
    dt = str2double( temp{2} );
    dt_all_h1(idNGA,1) = dt;
    gacc = fscanf(fid,'%15e ',npts);    
    gacc_all_h1{idNGA,1} = gacc;
    fclose(fid);
end

% Horizontal 2
npts_all_h2 = zeros(Ngm,1);
dt_all_h2 = zeros(Ngm,1);
gacc_all_h2 = cell(Ngm,1);
for idNGA = 1:Ngm
    clear gacc    
    currFn = at2fn_h2{idNGA}(2:end);
    fid = fopen(fullfile(RecordSetPath_orig,currFn),'r');
    for ll=1:3; fgets(fid); end; % Skip first 3 lines
    currLine = fgets(fid); % Get line containing data
    C = strsplit(currLine,', '); % Split line into two parts
    nptStr = C{1}; dtStr = C{2};
    temp = strsplit(nptStr,' '); % Further split string into parts
    npts = str2double( temp{2} );
    npts_all_h2(idNGA,1) = npts;
    temp = strsplit(dtStr,' '); % Further split string into parts
    dt = str2double( temp{2} );
    dt_all_h2(idNGA,1) = dt;
    gacc = fscanf(fid,'%15e ',npts);    
    gacc_all_h2{idNGA,1} = gacc;    
    fclose(fid);
end

% Vertical
npts_all_v = zeros(Ngm,1);
dt_all_v = zeros(Ngm,1);
gacc_all_v = cell(Ngm,1);
for idNGA = 1:Ngm
    clear gacc    
    currFn = at2fn_v{idNGA}(2:end);
    fid = fopen(fullfile(RecordSetPath_orig,currFn),'r');
    for ll=1:3; fgets(fid); end; % Skip first 3 lines
    currLine = fgets(fid); % Get line containing data
    C = strsplit(currLine,', '); % Split line into two parts
    nptStr = C{1}; dtStr = C{2};
    temp = strsplit(nptStr,' '); % Further split string into parts
    npts = str2double( temp{2} );
    npts_all_v(idNGA,1) = npts;
    temp = strsplit(dtStr,' '); % Further split string into parts
    dt = str2double( temp{2} );
    dt_all_v(idNGA,1) = dt;
    gacc = fscanf(fid,'%15e ',npts);    
    gacc_all_v{idNGA,1} = gacc;    
    fclose(fid);
end

%% Concisely store data for all components (for looping purposes)
newDirList = {'GMsX' 'GMsY' 'GMsZ'};
nComp = length(newDirList);
at2fn = [at2fn_h1 at2fn_h2 at2fn_v];
npts_all = [npts_all_h1 npts_all_h2 npts_all_v];
dt_all = [dt_all_h1 dt_all_h2 dt_all_v];
gacc_all = [gacc_all_h1 gacc_all_h2 gacc_all_v];

%% Save data extracted from NGAwest2 time series
save(fullfile(GMSMoutDir,timeSeriesOutMatFn),...
    'newDirList','at2fn','dt_all','npts_all','gacc_all');