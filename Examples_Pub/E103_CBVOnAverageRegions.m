% ------------------------------------------------------------------------- 
%                      E103_CBVOnLargeRegions
% 
% Performs deconvolution on growing regions on the synthetic flow model.
% Also displays results of the experiments.
% 
% 
% 
%                                          (c)Constantin Heck, 12-Mar-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 



clc;
close all;



%which results to show?
showCBVMaps = 1;
showPartition = 0;


%setup size of blocks to estimate CBV
blockSize = 10*[1,1];


%% load data

load smallDataSet.mat

%prepare sizes
m       = [size(Cmat,1),size(Cmat,2)];
n       = prod(m); %number of datapoints
mk      = [m,k]; %...

%do the downsampling
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



%% setup area where to run CBV-estimation

%initialize variables to store results
CBVest = zeros(m);
CBV    = zeros(m);
Cblock = zeros([m,klow]);

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

        
        %estimate CBV
        Cij   = Clow(idxi,idxj,:);
        Cij   = reshape(Cij,[],k);
        Cij   = mean(Cij,1)';
        CBVij = perfusion1c.cbvEstimation(Cij,timelinelow,AIFlow);
        
        
        %prepare Cbloc
        Cij = reshape(Cij,1,1,[]);
        Cij = repmat(Cij,cBlockSize(1),cBlockSize(2),1);

        

        %store results
        CBVest(idxi,idxj)   = CBVij;
        Cblock(idxi,idxj,:) = Cij;

        %fill the waitbar
        perc = (nBlocksj*(i-1) + j)/nBlocks;
        waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));       
        
        
        %prepare true flow to block-size
        CBVtrue           = phimat(idxi,idxj);
        CBVtrue           = mean(CBVtrue(:));
        CBV(idxi,idxj) = CBVtrue;
       
       
   end
end

delete(h);


%% showFlowMaps

if showCBVMaps
    
    %setup scaling cmax
    cScale = 1;
    
    %setup color axis
    tmp  = [CBV(:);CBVest(:)]*100;
    cmin = min(tmp);
    cmax = cScale*max(tmp);

    %setup 
    REMap    = abs(CBVest-CBV)./CBV;
    medianRE = median(REMap(:));

    %show flow
    figure(1);clf;
    colormap jet(512);
    set(gcf,'name','Results for CBV')


    subplot(1,3,1);
    imagesc(CBV*100);
    caxis([cmin,cmax]);
    axis image;
    title('True CBV (ml/100ml)')

    subplot(1,3,2);
    imagesc(CBVest*100);
    caxis([cmin,cmax]);
    axis image;
    title('Estimated CBV (ml/100ml)')

    %show division
    subplot(1,3,3);
    imagesc(REMap);
    caxis([0,1]);
    axis image;
    ti = sprintf('Relative Error in CBV |CBVest-CBVtrue|/CBVtrue\n medianRE=%1.4f',medianRE);
    title(ti)

end