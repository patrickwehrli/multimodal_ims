%% Intensity-based Image Registration 
% Includes manual control point (CPS) image registration, and automated
% intensity-based image registration.
%
% The script is to be run section wise to ensure correct data processing in
% following steps: (1) load reference images (e.g. PCA scores images)
%                  (2) Perform CPS registration for coarse alignment
%                  (3) SSIM-based optimization of registration parameters
%                  (4) Apply image registration on IMS dataset
%                  (5) save or load transformation matrix
%
%
%% clear all
clc;
clear all;
close all;
format long g;


%% Load reference images
[file,path] = uigetfile('*.*'); %select fixed image
cd(path)
fixed = imread(file);
fixed = mat2gray(fixed);
[frow, fcol] = size(fixed); %store size

[file,path] = uigetfile('*.*'); %select moving image
cd(path)
moving = imread(file);
moving = mat2gray(moving);


%% Manual control point registration
% for translational/rotational/shearing distorions; select at least 3 control points
cpselect(moving,fixed); % moving, fixed


%% Calculate the transformation matrix 
tformMat = fitgeotrans(movingPoints,fixedPoints,'affine'); % calculate transformation matrix
movingReg = imwarp(moving,tformMat,'OutputView',imref2d(size(moving)),'interp', 'bicubic'); %register moving image

initTform = tformMat; % define as initial for automatic registration

% Determine coarse alignment parameters 
fixedRefObj = imref2d(size(fixed)); % Default spatial referencing objects
movingRefObj = imref2d(size(moving));
[frow, fcol] = size(fixed);

% Apply coarse alignment
coarseReg = imwarp(moving, movingRefObj, initTform, 'OutputView', fixedRefObj, 'SmoothEdges', true);

% Display collage of overlays
figure
ax1 = subplot(2,2,1); imshow(fixed, []); title('fixed');
ax2 = subplot(2,2,2); imshow(coarseReg, []); title('moving after transform');
ax3 = subplot(2,2,3); imshowpair(fixed, coarseReg,'diff'); title('difference overlay');
ax4 = subplot(2,2,4); imshowpair(fixed, coarseReg,'ColorChannels','red-cyan'); title('red-cyan overlay');
linkaxes([ax1,ax2,ax3,ax4],'xy');


%% Display SSIM map
[ssimval2, ssimmap2] = ssim(fixed,coarseReg);
figure
imshow(ssimmap2,[]);
line4 = sprintf('(SSIM = %0.4f)',ssimval2);
title(line4);
truesize([frow*3 fcol*3]);
suptitle('manual coarse alignment')


%% SSIM-based optimization of registration parameters
% iteration through parameter settings
ssimMat = []; %preallocation
params = []; %preallocation

fixedRefObj = imref2d(size(fixed));
movingRefObj = imref2d(size(moving));
[frow, fcol] = size(fixed);

% set up registration metric
metric = registration.metric.MattesMutualInformation;
metric.NumberOfSpatialSamples = 500; %default is 500
metric.NumberOfHistogramBins = 60; %default is 50
metric.UseAllPixels = true; %default is true

optimizer = registration.optimizer.OnePlusOneEvolutionary;
optimizer.MaximumIterations = 100; %default is 100


for gf = 1.005:0.0005:1.1
    optimizer.GrowthFactor = gf;
   
    for oe = 1.5e-07:1:1.5e-03
        optimizer.Epsilon = oe;
        
        for ir = 6.25e-04:1:6.25e-02
            optimizer.InitialRadius = ir;
        
            % Apply transformation
            tformMat = imregtform(moving,movingRefObj,fixed,fixedRefObj,'affine',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
            %tformMat = imregtform(moving,movingRefObj,fixed,fixedRefObj,'affine',optimizer,metric,'PyramidLevels',3);
            movingReg = imwarp(moving, movingRefObj, tformMat, 'OutputView', fixedRefObj,'interp', 'bilinear', 'SmoothEdges', true);

            % SSIM 
            [ssimVal, ssimMap] = ssim(fixed,movingReg);
            ssimMat = [ssimMat ssimVal];
            params = [params; gf oe ir];
        
        end
    end
end

[maxSSIM idx] = max(ssimMat); % parameters with maxi SSIM value

optimizer.GrowthFactor = params(idx,1);
optimizer.Epsilon = params(idx,2);
optimizer.InitialRadius = params(idx,3);

% run registration with new parameters
tformMat = imregtform(moving,movingRefObj,fixed,fixedRefObj,'affine',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
movingReg = imwarp(moving, movingRefObj, tformMat, 'OutputView', fixedRefObj,'interp', 'bilinear', 'SmoothEdges', true);

% Display SSIM map
[ssimval2, ssimmap2] = ssim(fixed,movingReg);
sprintf('(SSIM = %0.4f)',ssimval2)

% Jaccard index calculation
im1 = fixed; im1 = im2uint8(im1); im1=double(im1)+1;im1(im1>255)=255;
im2 = movingReg; im2 = im2uint8(im2); im2=double(im2)+1;im2(im2>255)=255;
JaccardIdx1 = jaccard(im1,im2); JaccardIdx1(isnan(JaccardIdx1))=0;
JaccardSum = sum(JaccardIdx1);
sprintf('(Jaccard index = %0.4f)',JaccardSum)

% Mutual information calculation
im1 = fixed;im1 = im2uint8(im1); [~,~,indrow] = unique(im1(:));
im2 = movingReg;im2 = im2uint8(im2);[~,~,indcol] = unique(im2(:));
jointHistogram = accumarray([indrow indcol], 1);
jointProb = jointHistogram / numel(indrow);
indNoZero = jointHistogram ~= 0;
jointProb1DNoZero = jointProb(indNoZero);
jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));
mutualInformation = entropy1 + entropy2 - jointEntropy;
sprintf('(Mutual information = %0.4f)',mutualInformation)


%% Display collage of overlays
figure
ax1 = subplot(2,2,1); imshow(fixed, []); title('fixed');
ax2 = subplot(2,2,2); imshow(movingReg, []); title('moving after transform');
ax3 = subplot(2,2,3); imshowpair(fixed, movingReg,'diff'); title('difference overlay');
ax4 = subplot(2,2,4); imshowpair(fixed, movingReg,'ColorChannels','red-cyan'); title('red-cyan overlay');
linkaxes([ax1,ax2,ax3,ax4],'xy');


%% Display SSIM map
[ssimval2, ssimmap2] = ssim(fixed,movingReg);
figure
imshow(ssimmap2,[]);
line4 = sprintf('(SSIM = %0.4f)',ssimval2);
title(line4);
truesize([frow*3 fcol*3]);
suptitle('final registration')


%% apply transformation on 3D IMS dataset
dset = dset; % provide appropriate dataset
dsetReg = imwarp(dset, tformMat, 'nearest', 'OutputView', imref2d(size(fixed)), 'SmoothEdges', true);
dsetReg = double(dsetReg);
figure; imshow(dsetReg(:,:,1),[])


%% save transformation matrix to file
[filename, pathname, filterindex] = uiputfile( ...
    {'*.mat','MAT-files (*.mat)'},...
    'Save file as',...
    'tformMat');
cd(pathname);
save(filename,'tformMat');      


%% load transformation matrix
[FileName,PathName] = uigetfile('*.mat','Select file');
cd(PathName);
Pstr = load(FileName);
tformMat = getfield(Pstr, 'tformMat');



