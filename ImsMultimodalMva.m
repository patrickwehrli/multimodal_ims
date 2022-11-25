%% Multimodal IMS Multivariate Analysis Script
% This script concatenates (registered) IMS matrices and exports them for 
% multivariate data analysis in SIMCA software. In a second step scores
% images from MVA are reconstructed from SIMCA scores export.
% 
%  (1a) Supply at least one dataset and m/z vector for image MVA
%  (1b) For more than one dataset, the matrices will be concatenated
%  (1c) Black pixels i.e., pixels from e.g. off-tissue IMS analysis will be
%       removed and stored (to be later recombined)
%  (1d) Datasets to be analyzed in SIMCA software are saved to file.
%  (2a) Scores from SIMCA are imported and recombined with black pixels
%  (2b) Scores images are reconstructed and saved to file with and without
%       titles

% allocate datasets to be combine
dset1 = [dsetPeakData]; % data matrix modality 1
dset2 = [];             % data matrix modality 2
dset3 = [];             % data matrix modality 3

mzs1 = pksMzsFullRes.'; % mz vector modality 1
mzs2 = [];              % mz vector modality 2
mzs3 = [];              % mz vector modality 3

nPixels = size(dset1,1); % create pixel IDs
ID = transpose(1:nPixels);
dsetCombined = [ID, dset1, dset2, dset3]; % concatenate ID and datasets

% Separate data pixels from pixels without data (i.e. black pixels) 
dsetCombined(isnan(dsetCombined))=0;
dsetCombined(dsetCombined<=0)=0;
S = sum(dsetCombined(:,2:end),2);
idxData = S(:,1) > 0;
idxBlack = S(:,1) <= 0;

dataPixels = dsetCombined(idxData,:);
blackPixels = dsetCombined(idxBlack,:); 

% create headers
mzsMod1s = string(mzs1.'); mzsMod1s = strcat('mod1_', mzsMod1s);
mzsMod2s = string(mzs2.'); mzsMod2s = strcat('mod2_', mzsMod2s);
mzsMod3s = string(mzs3.'); mzsMod3s = strcat('mod3_', mzsMod3s);

headers = ['ID', mzsMod1s, mzsMod2s, mzsMod3s];

% check for duplicates in header
[~, uniqueIdx] = unique(headers); % indeces of unique headers
duplicates = headers; 
duplicates(uniqueIdx) = []; % delete unique ones to get duplicates only
duplicates = unique(duplicates); % find the unique duplicates
nDuplicates = length(duplicates); % number of duplicates

% warn for duplicates
disp(['number of duplicates in headers: [' num2str(nDuplicates) ']']) ;

% export matrices
[name, path] = uiputfile('combinedDataset.txt','Save as');
cd(path)

% write strings to file
fid = fopen(name,'wt');
    if fid ~= -1
        fprintf(fid,'%s\t',headers);
        fprintf(fid,'\r\n');
        fclose(fid);
    end

dlmwrite(name,dataPixels,'-append','delimiter','\t','precision',16,'newline','pc'); % append dataPixels to file
dlmwrite('blackPixels.txt',blackPixels,'delimiter','\t','precision',16,'newline','pc'); % write blackPixels to file

% display the sizes of the output matrices
disp(['blackPixels size: [' num2str(size(blackPixels)) ']']) ;
disp(['dataPixels size: [' num2str(size(dataPixels)) ']']) ;
disp(['pixel sum: [' num2str(size(dataPixels,1)+size(blackPixels,1)) ']' ' should be [' num2str(size(dsetCombined,1)) ']' ]) ;


%% Score image reconstruction from SIMCA scores 
% SIMCA software scores saved as .txt file in SIMCA are reconstructed to
% scores images.

scores = readtable('scores.txt'); % read text file
scores = table2array(scores);
imgHeight = 100; % edit image height

% Recombine scores and black pixels
nScores = size(scores,2);
if isempty(blackPixels) == false
    blackPixels = blackPixels(:,1:nScores);
    blackPixels = double(blackPixels);
    scoresComplete = vertcat(scores, blackPixels);
else
    scoresComplete = scores;
end

scoresComplete = sortrows(scoresComplete,1); % sort combined matrix after ID
score1 = reshape(scoresComplete(:,2),imgHeight,[]);score1 = mat2gray(score1); % preview first scores image

cmap = hot(256); %colormap
figure; imshow(score1,'Colormap',cmap,'InitialMagnification','fit'); title('preview');
colorbar

% Save all PCA scores images to file
scoresDir = uigetdir('C:\');
cd(scoresDir);

% Save images with titles
for i = 2:nScores
    t = i-1;
    iScore = reshape(scoresComplete(:,i),imgHeight,[]); iScore = mat2gray(iScore);
    iImg = figure('visible','off');
    imshow(iScore,'Colormap',cmap,'InitialMagnification','fit'); 
    iTitle = strcat('scores ',num2str(t));
    title(iTitle);
    
    outputName = sprintf('scores_%d.png', t);
    saveas(iImg,outputName,'png')
end

% Save images without margins (frame and title)
if exist('scores images without margins','dir') == 0
    mkdir('scores images without margins');
end

cd('scores images without margins')

for i = 2:nScores
    t = i-1;
    iScore = reshape(scoresComplete(:,i),imgHeight,[]);
    iScore = mat2gray(iScore);
    IndScore = gray2ind(iScore,256);
    
    outputName = sprintf('scores_t%d.png', t);
    imwrite(IndScore,cmap,outputName);
 end
