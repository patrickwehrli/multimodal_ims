%% MS peak extraction routine  
% This script performs peak picking and peak data extraction at full
% mass resolution. The purpose is to reduce data size while retaining the
% number of variables to obtain a dataset suitable for multivariate
% analysis.
% 
% The script is to be run section wise to ensure correct data processing in
% following steps: (1) Generate baseline corrected sum spectrum, 
%                  (2) Automated peak picking
%                  (3) Manually adjust bin limits
%                  (4) save bin limits to file (if desired)
%                  (5) read bin limits from file (if desired)
%                  (6) Peak extraction at full mass resolution
%                  (7) Check for duplicate variables
%                  (8) Remove duplicate variables
%

%% Generate baseline corrected sum spectrum

dsetX = dset(:,2:end).'; % provide MS dataset [mz x pixels]
mzs = dset(:,1); % extracts mz vector

sumSpect = sum(dsetX).'; % sum spectrum and transpose 

% baseline correction; Edit parameters to suit data!
sumSpect_baselined = msbackadj(mzs, sumSpect, 'WindowSize',30,'StepSize',...
    30,'RegressionMethod','pchip','EstimationMethod','quantile',...
    'SmoothMethod','none','Quantile',0.01, 'SHOWPLOT',true);
    
figure; plot(mzs,sumSpect_baselined); % plot baseline-corrected sum spectrum


%% Automated peak picking
% peak maxima, location indices of peaks (if mz's are provided locs are mz values)

[pks,locs,w,p] = findpeaks(sumSpect_baselined,'MinPeakProminence',4); % find peaks

% determine bin start positions for peak extraction
binStartIdx = locs-w; 
binStartIdx(binStartIdx<0)=1;
binStartIdx = round(binStartIdx);
binStartmzs = mzs(binStartIdx);

% determine bin stop positions for peak extraction
binStopIdx = locs+1.6*w; 
binStopIdx = round(binStopIdx);
binStopmzs = mzs(binStopIdx);

pksMzs =  mzs(locs); % get mz values of peaks
figure; plot(mzs,sumSpect_baselined,pksMzs,pks,'kv','markerfacecolor',[0 0 0], 'MarkerSize', 3); % plot spectrum

% plot bin limits
for k = 1:length(binStartmzs)
    line([binStartmzs(k,:) binStartmzs(k,:)],[-100 p(k,:)]);
    line([binStopmzs(k,:) binStopmzs(k,:)],[-100 p(k,:)]);
end

binLimitIdx = [binStartIdx binStopIdx]; % Get start and stop indeces of bins

fprintf('%.0f peaks picked',size(pks,1)) % print number of peaks picked


%% Manually adjust bin limits
% Apply if needed
binLimits = [mzs(binStartIdx).' mzs(binStopIdx).']; % start and stop m/z values of bins (row)
msviewer(mzs,sumSpect_baselined,'Markers', binLimits); % MS-viewer window

f = warndlg({'Adjust bin limits and export as "x",','after that click ok!'},'Bin Limits');
waitfor(f);

[~,binLimitIdx] = min(abs(mzs-x.')); % get closest value to selected peak the mz matrix
binLimitIdx = reshape(binLimitIdx,2,[]).'; % new start and stop indeces of bins

NewBinLimitIdx = reshape(binLimitIdx,1,[]); % reshape for plotting
NewBinLimits = mzs(NewBinLimitIdx); % reshape for plotting
msviewer(mzs,sumSpect_baselined,'Markers', NewBinLimits); % plot for control


%% save bin limits to file (if desired)
Ydir = uigetdir('C:\');
cd(Ydir);

dlmwrite('BinLimitIdx.txt',binLimitIdx,'precision',7); % save bin limit indeces to file (number of significant digits)
dlmwrite('BinPeaksMzs.txt',pksMzs,'precision',7);% save peak's top m/z values to file (required for binned data)
mzLims =  mzs(binLimitIdx); 
dlmwrite('BinLimitMz.txt',mzLims,'precision',7); % save bin limits as m/z values


%% read bin limits from file (if desired)
[FileName,FilePath] = uigetfile('*.txt'); % browse to file location
cd(FilePath);

binLimitIdx = readtable(FileName,'ReadVariableNames', false); %parse bin data
binLimitIdx = table2array(binLimitIdx);

% display bins
NewBinLimitIdx = reshape(binLimitIdx,1,[]);
NewBinLimits = mzs(NewBinLimitIdx);
msviewer(mzs,sumSpect,'Markers', NewBinLimits)

%% Peak extraction at full mass resolution 
% extracts peak data within bin limits in full mass resolution
IdxVec = [1:size(dset4bin,2)]; % mass range index vector
nBins = size(binLimitIdx,1); % number of bins
nPix = size(dset4bin,1); % number of pixels

% count number of total data points
nDataPoints = 0;
for j = 1:nBins 
    n = sum(IdxVec >= binLimitIdx(j,1) & IdxVec <= binLimitIdx(j,2));
    nDataPoints = nDataPoints + n;
end

% collect data points within bin limits
dsetPeakData = zeros(nPix,nDataPoints); % preallocation
for i = 1:nPix 
    peakData = [];
    for j = 1:nBins 
        jBin = dset4bin(i, find(IdxVec >= binLimitIdx(j,1) & IdxVec <= binLimitIdx(j,2))); % extract values of each bin
        peakData = [peakData jBin]; % combine peak values for each pixel
    end
    dsetPeakData(i,:) = peakData; % add all peak values of each pixel to a combined matrix
end

% create a corresponding m/z vector
pksMzsFullRes =[];
for j = 1:nBins
    jMz = mzs(find(IdxVec >= binLimitIdx(j,1) & IdxVec <= binLimitIdx(j,2)),1); % extract m/z values
    pksMzsFullRes = [pksMzsFullRes; jMz]; % combine m/z values of each datapoint
end

size(dsetPeakData) % display array size

% display bar graph of peak data at full mass resolution
figure
bar(pksMzsFullRes,dsetPeakData(1000,:)); %If no peaks are visible, zoom in or choose different pixel! 


%% Check for duplicate variables
% Overlapping bin limits create duplicate variables which need to be
% removed for multivariate analysis
[~, uniqueIdx] = unique(pksMzsFullRes); % indeces of unique headers
duplicates = pksMzsFullRes; 
duplicates(uniqueIdx) = []; % delete unique ones to get duplicates only
duplicates = unique(duplicates); % find the unique duplicates
nDuplicates = length(duplicates) % display number of duplicates


%% Remove duplicate variables
headers = pksMzsFullRes;
dsetX = dsetPeakData;

szDsetX = size(dsetX,2); % dset size
[uniqueHeaders, uniqueIdx] = unique(headers); % get unique headers and their indeces
dsetX = dsetX(:,uniqueIdx); % extract unique variables

% display the sizes of the output matrices
disp(['headers size: [' num2str(size(headers,1)) ']']);
disp(['unique headers size: [' num2str(size(uniqueHeaders,1)) ']' ' should be [' num2str(size(headers,1)-nDuplicates) ']' ]) ;

disp(['dset size: [' num2str(szDsetX) ']']) ;
disp(['unique dset size: [' num2str(size(dsetX,2)) ']' ' should be [' num2str(size(uniqueHeaders,1)) ']' ]) ;

pksMzsFullRes = uniqueHeaders; % rename output
dsetPeakData = dsetX; % rename output


