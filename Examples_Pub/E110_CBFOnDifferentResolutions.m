function [] = E110_CBFOnDifferentResolutions
% ------------------------------------------------------------------------- 
%                      E110_CBFOnDifferentResolutions
% 
% Performs deconvolution on different resolutions on the synthetic flow 
% model.
% Also displays results of the experiments.
% 
% 
% 
%                                          (c)Constantin Heck, 28-Feb-2017 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clc;
close all;
clear;


%setup main parameters
OIList = 0.01*[1, 1, 1, 1, 1, 1, 1];
voxelSizeFactorList = [1,2,4,8,16,32,64];


% load data

D = load('smallDataSet.mat');
Cmat = D.Cmat;
Fmat = D.Fmat;
qmat = D.qmat;
aifval = D.aifval;
k = D.k;
perfmat = D.perfmat;
timeline = D.timeline;
prm = D.prm;

%prepare sizes
m       = [size(Cmat,1),size(Cmat,2)];
% n       = prod(m); %number of datapoints
% mk      = [m,k]; %...

%do the downsampling
step        = 1;
idxsamp     = (1:step:k);
Clow        = squeeze(Cmat(:,:,idxsamp));
timelinelow = timeline(idxsamp);
AIFlow      = aifval(idxsamp);

klow    = numel(idxsamp);
% mklow   = [m,klow];
nVoxelSizes = length(voxelSizeFactorList);

%store relative errors
zeroVec = zeros(nVoxelSizes, 1);
% meanREs = struct('DeconvPerf', zeroVec, 'MSPerf', zeroVec, 'DeconvPerfLoc', zeroVec, 'MSPerfLoc', zeroVec);
% stdREs = struct('DeconvPerf', zeroVec, 'MSPerf', zeroVec, 'DeconvPerfLoc', zeroVec, 'MSPerfLoc', zeroVec);

% Run for different voxel sizes

for vs = (1:nVoxelSizes)
    
    %setup an empty image do display the blocks
    voxelSizeFactor = voxelSizeFactorList(vs);
    OI = OIList(vs);

    %determine number of blocks in i and j direction, determine blocksize
    blockSize = voxelSizeFactor*[1,1];
    partitionImage = zeros(m(1),m(2));
    nBlocksi = ceil(m(1)/blockSize(1));
    nBlocksj = ceil(m(2)/blockSize(2));
    nBlocks  = nBlocksi*nBlocksj;
%     numBlock = prod(blockSize);

    % setup area where to run the deconvolution
    %initialize variables to store resultss
    CBFCirc = zeros(m);
    CBFMS = zeros(m);
    Iest = zeros([m,2*klow]);
    Cest = zeros([m,2*klow]);
    Cblock = zeros([m,klow]);
    Fsource = Fmat; Fsource(Fsource<0)=0;
    P = zeros(m);
    PLocal = zeros(m);

    % get deconvolution matrix
    deltaT  = timelinelow(2)-timelinelow(1);
    A       = perfusion1c.getCircularConvolutionMatrix(AIFlow,deltaT);
    fprintf('Starting SVD...');
    tic; [U,S,V] = svd(A);
    fprintf('...done. Elapsed time: %1.3fs\n',toc);

    %do the job in blocks
    h = waitbar(0);
    for i = 1:nBlocksi
       msg = ['Block x: ' int2str(i) ', ' int2str(nBlocksi)];
       disp(msg);
       for j = 1:nBlocksj
           
            %indices of current block
            idxi = blockSize(1)*(i-1) + (1:blockSize(1));
            idxj = blockSize(2)*(j-1) + (1:blockSize(2));
            %remove indices larger then maximal entries
            idxi(idxi>m(1)) = [];
            idxj(idxj>m(2)) = [];
            %get actual size of current block
            cBlockSize = [numel(idxi),numel(idxj)];
%             cBlockn    = prod(cBlockSize);
            %create image displaying the block
            partitionImage(idxi,idxj) = (i-1)*nBlocksj + j;
            %prepare data curve
            Cij = Clow(idxi,idxj,:);
            Cij = reshape(Cij,[],k);
            Cij = mean(Cij,1)';
            %estimate flow
            [FCirc,Irec,Crec] = perfusion1c.circularDeconvolution(Cij,timelinelow,OI,U,S,V);
            FMS               = perfusion1c.maximumSlope(Cij,timelinelow,AIFlow);
            %prepare Irec and Crec to store them
            Irec = reshape(Irec,1,1,[]);
            Irec = repmat(Irec,cBlockSize(1),cBlockSize(2),1);
            Crec = reshape(Crec,1,1,[]);
            Crec = repmat(Crec,cBlockSize(1),cBlockSize(2),1);
            %prepare Cblock, averaged concentrations
            Cij = reshape(Cij,1,1,[]);
            Cij = repmat(Cij,cBlockSize(1),cBlockSize(2),1);
            %store results
            CBFCirc(idxi,idxj) = FCirc;
            CBFMS(idxi,idxj) = FMS;
            Iest(idxi,idxj,:) = Irec;
            Cest(idxi,idxj,:) = Crec;
            Cblock(idxi,idxj,:) = Cij;
            %fill the waitbar
            perc = (nBlocksj*(i-1) + j)/nBlocks;
            waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));       
            %prepare true perfusion for block
            Ftrue = perfmat(idxi,idxj); %perfmat is the streamline perfusion
            Ftrue = mean(Ftrue(:));
            P(idxi,idxj) = Ftrue;
            %prepare true local perfusion for block
            Ftop = qmat{1}(idxi(1),idxj);
            Fleft = qmat{2}(idxi,idxj(1));
            Fbdry = sum(Ftop)+sum(Fleft);
            Fs = sum(reshape(Fsource(idxi,idxj),[],1));
            hdi = numel(idxi)*numel(idxj)*prod(prm.h);
            PLocal(idxi,idxj) = (Fbdry + Fs)/hdi; %total inflow divided by volume
       end
    end
    delete(h);
    
    % calculate relative errors and store average perfusion
    
    %setup relative errors with respect to streamline perfusion
    results.CBFCirc{vs} = CBFCirc;
    results.CBFMS{vs} = CBFMS;
    results.P{vs} = P;
    results.PLocal{vs} = PLocal;
%     RECirc = abs(CBFCirc-P)./P;
%     REMS = (CBFMS-P)./P;
%     results.DeconvPerf{vs} = RECirc(:);
%     results.MSPerf{vs} = REMS(:);
%     stdREs.DeconvPerf = std(RECirc(:));
%     stdREs.MSPerf = std(REMS(:));
    
    %setup relative errors with respect to local perfusion
%     RECirc = (CBFCirc-PLocal)./PLocal;
%     REMS = (CBFMS-PLocal)./PLocal;    
%     results.DeconvPerfLoc{vs} = RECirc(:);
%     results.MSPerfLoc{vs} = REMS(:);
%     stdREs.DeconvPerfLoc = std(RECirc(:));
%     stdREs.MSPerfLoc = std(REMS(:));
    
end

results.voxelSizeFactorList = voxelSizeFactorList;
results.prm = prm;
pathsave = [mfilename '-' 'results.mat'];
save(pathsave,'results');
msg = ['Saving ' pathsave];
disp(msg);
% % show results goes into a separate file
% 
% figure(1)
% suplot(2,2,1)
% plot(voxelSizeFactorList, meanREs.DeconvPerf)
% title('meanRE: Deconv-Perf')
% suplot(2,2,2)
% plot(voxelSizeFactorList, meanREs.MSPerf)
% title('meanRE: MS-Perf')
% suplot(2,2,3)
% plot(voxelSizeFactorList, meanREs.DeconvPerfLoc)
% title('meanRE: Deconv-PerfLoc')
% suplot(2,2,4)
% plot(voxelSizeFactorList, meanREs.MSPerfLoc)
% title('meanRE: MS-PerfLoc')

% figure(1)
% suplot(2,2,1)
% plot(voxelSizeFactorList, stdREs.DeconvPerf)
% title('meanRE: Deconv-Perf')
% suplot(2,2,2)
% plot(voxelSizeFactorList, stdREs.MSPerf)
% title('meanRE: MS-Perf')
% suplot(2,2,3)
% plot(voxelSizeFactorList, stdREs.DeconvPerfLoc)
% title('meanRE: Deconv-PerfLoc')
% suplot(2,2,4)
% plot(voxelSizeFactorList, stdREs.MSPerfLoc)
% title('meanRE: MS-PerfLoc')