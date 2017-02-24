% ------------------------------------------------------------------------- 
%                      E102_CBFOnLargeRegions
% 
% Performs deconvolution on growing regions on the synthetic flow model.
% Also displays results of the experiments.
% 
% 
% 
%                                          (c)Constantin Heck, 24-Feb-2017 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clc;
close all;
clear;

%which results to show?
compareWithStreamlinePerfusion = 1;
compareWithLocalPerfusion = 1;
showSingleCurve = 0;

%setup main parameters
blockSize = 10*[1,1];
OI        = .01; %oscillation index (threshold for bSVD reconstruction)


%% load data

load smallDataSet.mat

%prepare sizes
m       = [size(Cmat,1),size(Cmat,2)];
n       = prod(m); %number of datapoints
mk      = [m,k]; %...

%do the downsampling
step        = 1;
idxsamp     = (1:step:k);
Clow        = squeeze(Cmat(:,:,idxsamp));
timelinelow = timeline(idxsamp);
AIFlow      = aifval(idxsamp);

klow    = numel(idxsamp);
mklow   = [m,klow];



%% prepare blocks

%setup an empty image do display the blocks
partitionImage = zeros(m(1),m(2));

%determine number of blocks in i and j direction
nBlocksi = ceil(m(1)/blockSize(1));
nBlocksj = ceil(m(2)/blockSize(2));
nBlocks  = nBlocksi*nBlocksj;
numBlock = prod(blockSize);


%% setup area where to run the deconvolution


%initialize variables to store results
CBFCirc = zeros(m);
CBFMS   = zeros(m);
Iest    = zeros([m,2*klow]);
Cest    = zeros([m,2*klow]);
Cblock  = zeros([m,klow]);
Fsource = Fmat; Fsource(Fsource<0)=0;

%get longer timeline
timelinelow2 = linspace(0,2*timelinelow(end),2*klow)';

%get deconvolution matrix
deltaT  = timelinelow(2)-timelinelow(1);
A       = perfusion1c.getCircularConvolutionMatrix(AIFlow,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);

%do the job in blocks
h = waitbar(0);
for i = 1:nBlocksi
   for j = 1:nBlocksj
       
        %setup current block
        idxi = blockSize(1)*(i-1) + (1:blockSize(1));
        idxj = blockSize(2)*(j-1) + (1:blockSize(2));
       
        %remove indices larger then maximal entries
        idxi(idxi>m(1)) = [];
        idxj(idxj>m(2)) = [];
        
        %get current blockSize
        cBlockSize = [numel(idxi),numel(idxj)];
        cBlockn    = prod(cBlockSize);
       
        
        %create image to display the blocks
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
        
        %prepare Cbloc
        Cij = reshape(Cij,1,1,[]);
        Cij = repmat(Cij,cBlockSize(1),cBlockSize(2),1);
        
        %store results
        CBFCirc(idxi,idxj)  = FCirc;
        CBFMS(idxi,idxj)    = FMS;
        Iest(idxi,idxj,:)   = Irec;
        Cest(idxi,idxj,:)   = Crec;
        Cblock(idxi,idxj,:) = Cij;

        %fill the waitbar
        perc = (nBlocksj*(i-1) + j)/nBlocks;
        waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));       
        
        
        %prepare true perfusion for block
        Ftrue          = perfmat(idxi,idxj);
        Ftrue          = mean(Ftrue(:));
        P(idxi,idxj)   = Ftrue;
        
        %prepare true local perfusion for block
        Ftop              = qmat{1}(idxi(1),idxj); 
        Fleft             = qmat{2}(idxi,idxj(1)); 
        Fbdry             = sum(Ftop)+sum(Fleft);
        Fs                = sum(reshape(Fsource(idxi,idxj),[],1));
        hdi               = numel(idxi)*numel(idxj)*prod(prm.h);
        PLocal(idxi,idxj) = (Fbdry + Fs)/hdi;
       
       
   end
end

delete(h);


%% Relative Error vs Perfusion

if compareWithStreamlinePerfusion
    
    %setup scaling for plotting with cmax
    cScale   = 1;
    
    %setup true perfusion and relative errors
    RECirc  = (CBFCirc-P)./P;
    REMS    = (CBFMS-P)./P;    
    mRECirc = mean(RECirc(:));
    mREMS   = mean(REMS(:));
    sRECirc = std(RECirc(:));
    sREMS   = std(REMS(:));
    
    %setup color axis
    tmp    = [P(:);CBFCirc(:)]*100*60;
    caCirc = [min(tmp),cScale*max(tmp)];
    tmp    = [P(:);CBFMS(:)]*100*60;
    caMS   = [min(tmp),cScale*max(tmp)];   

    %results
    figure(1);clf;
    colormap parula;
    set(1,'name','PERFUSION');

    
    %CBF and CBFCirc
    subplot(2,3,1);
    imagesc(P*100*60);
    caxis(caCirc);
    axis image;
    title('Perfusion (ml/min/100ml)')
    
    subplot(2,3,2);
    imagesc(CBFCirc*100*60);
    caxis(caCirc);
    axis image;
    title('Estimated flow Circ (ml/min/100ml)')

    %show division
    subplot(2,3,3);
    imagesc(RECirc);
    axis image;
    ti = sprintf('Relative Errors Circ\nmeanRE: %1.4f, stdRE=%1.4f\n',mRECirc,sRECirc);
    title(ti)
    
    %CBF MS
    subplot(2,3,4);
    imagesc(P*100*60);
    caxis(caMS);
    axis image;
    title('Perfusion (ml/min/100ml)')
    
    subplot(2,3,5);
    imagesc(CBFMS*100*60);
    caxis(caMS);
    axis image;
    title('Estimated flow MS (ml/min/100ml)')

    %show division
    subplot(2,3,6);
    imagesc(REMS);
    axis image;
    ti = sprintf('Relative Errors MS\nmeanRE: %1.4f, stdRE=%1.4f\n',mREMS,sREMS);
    title(ti)    
    

end





%% relative error vs local perfusion


if compareWithLocalPerfusion
    
    %setup scaling cmax
    cScale   = 1;
    
    %setup true perfusion and relative errors
    RECirc   = (CBFCirc-PLocal)./PLocal;
    REMS     = (CBFMS-PLocal)./PLocal;    
    mRECirc  = mean(RECirc(:));
    sRECirc  = std(RECirc(:));
    mREMS    = mean(REMS(:));
    sREMS    = std(REMS(:));
    
    %setup color axis
    tmp    = [PLocal(:);CBFCirc(:)]*100*60;
    caCirc = [min(tmp),cScale*max(tmp)];
    tmp    = [PLocal(:);CBFMS(:)]*100*60;
    caMS   = [min(tmp),cScale*max(tmp)];   

    %results
    figure(2);clf;
    colormap parula;
    set(2,'name','LOCAL PERFUSION');

    
    %CBF Circ
    subplot(2,3,1);
    imagesc(PLocal*100*60);
    caxis(caCirc);
    axis image;
    title('Local Perfusion (ml/min/100ml)')
    
    subplot(2,3,2);
    imagesc(CBFCirc*100*60);
    caxis(caCirc);
    axis image;
    title('Estimated flow Circ (ml/min/100ml)')

    %show division
    subplot(2,3,3);
    imagesc(RECirc);
    axis image;
    ti = sprintf('Relative Errors Circ\nmeanRE: %1.4f, stdRE=%1.4f\n',mRECirc,sRECirc);
    title(ti)
    
    %CBF MS
    subplot(2,3,4);
    imagesc(PLocal*100*60);
    caxis(caMS);
    axis image;
    title('Local Perfusion (ml/min/100ml)')
    
    subplot(2,3,5);
    imagesc(CBFMS*100*60);
    caxis(caMS);
    axis image;
    title('Estimated flow MS (ml/min/100ml)')

    %show division
    subplot(2,3,6);
    imagesc(REMS);
    axis image;
    ti = sprintf('Relative Errors MS\nmeanRE: %1.4f, stdRE=%1.4f\n',mREMS,sREMS);
    title(ti)    
    

end





%% results for a single voxel

if showSingleCurve
    
    %setup squeeze function
    s = @(v) squeeze(v(:));    
    
    %position of single curve
    pos = [32,30];
%     pos = randi(m(1),2,1);
    
    %get the stuff
    corig = s(Cblock(pos(1),pos(2),:));
    cest  = s(Cest(pos(1),pos(2),:));
    iest  = s(s(Iest(pos(1),pos(2),:)));
    
    %double-check
    tau   = timelinelow(2)-timelinelow(1);
    ctest = conv(AIFlow,iest);
    ctest = tau*ctest(1:k);


    figure(3);clf;
    set(3,'name','Comparison on single voxel');

    subplot(1,3,1);
    colormap gray(512);
    imagesc(squeeze(Cblock(:,:,20)))
    axis image;
    hold on;
    plot(pos(1),pos(2),'xb')
    hold off;
    title(sprintf('Ctrue (blue) and reconstructed C (red) at position pos=%i,%i',pos(1),pos(2)))    
    
    subplot(1,3,2);
    plot(timelinelow,corig,timelinelow,ctest,'LineWidth',3);
    title(sprintf('Ctrue (blue) and reconstructed C (red) at position pos=%i,%i',pos(1),pos(2)))
    
    subplot(1,3,3);
    plot(timelinelow,iest(1:k),'LineWidth',3);
    title('estimated I at position pos')
    

end