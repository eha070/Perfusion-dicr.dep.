% ------------------------------------------------------------------------- 
%                      E06_DeconvolutionOnGrowingRegions 
% 
% Performs deconvolution on growing regions on the synthetic flow model.
% 
% 
% 
%                                          (c)Constantin Heck, 12-Mar-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 



clc;
close all;



%setup which flow-calculation to use
indicatorcalc = 'conv';
% indicatorcalc = 'PDE';



%setup the flow to compare the results to
trueFlow = 'perfusion';
% trueFlow = 'averaging'; %do NOT use this for indicatorcalc='conv'



%which results to show?
showCBVMaps        = 1;
writeData          = 1;
showPartition      = 0;


%setup size of blocks to run the deconvolution on
blockSize = 10*[1,1];


%prepare downsampling of data (shorter runtime for SVD)
%It will hold: Clow = Clow(:,:,1,1:step:k);
step = 1;





%% load data

%settings
prm               = settings;
basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
basenameFlow      = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');
foldername        = './results/';
csi               = prm.csi;



%setup paths
pathloadFlow = ['synt-createflowTPFA-' basenameFlow '.mat'];
switch indicatorcalc
    case 'PDE'
        pathload = [foldername,'synt-createindicatorpde-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
    case 'conv'
        pathload = [foldername,'syntconv-createindicatorconv-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
    otherwise
        error('Method set for indicatorcalc unknown')
end

% load concentration curves
msg = ['Loading ' pathload];
disp(msg);
D = load(pathload);

%setup C etc
Cmat     = squeeze(D.Cmat);
AIF      = D.aifval;
timeline = D.prm.timeline;

%setup CBF
msg = ['Loading ' pathloadFlow];
disp(msg);
E = load(pathloadFlow);

CBV = E.phimat;
%clear memory
clearvars E D;




%% prepare data
%prepare sizes
m       = [size(Cmat,1),size(Cmat,2)];
n       = prod(m); %number of datapoints
k       = size(Cmat,3); %number of timepoints
mk      = [m,k]; %...

%do the downsampling
idxsamp     = (1:step:k);
Clow        = squeeze(Cmat(:,:,idxsamp));
timelinelow = timeline(idxsamp);
AIFlow      = AIF(idxsamp);

klow    = numel(idxsamp);
mklow   = [m,klow];



%% prepare blocks

%setup an empty image do display the blocks
partitionImage = zeros(csi(1),csi(2));

%determine number of blocks in i and j direction
nBlocksi = ceil(csi(1)/blockSize(1));
nBlocksj = ceil(csi(2)/blockSize(2));
nBlocks  = nBlocksi*nBlocksj;
numBlock = prod(blockSize);


%% setup area where to run CBV-estimation



%initialize variables to store results
CBVest = zeros(m);
Cblock = zeros([m,klow]);

%do the job in blocks
h = waitbar(0);
for i = 1:nBlocksi
   for j = 1:nBlocksj
       
        %setup current block
        idxi = blockSize(1)*(i-1) + (1:blockSize(1));
        idxj = blockSize(2)*(j-1) + (1:blockSize(2));
       
        %remove indices larger then maximal entries
        idxi(idxi>csi(1)) = [];
        idxj(idxj>csi(2)) = [];
        
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
        CBVtrue        = CBV(idxi,idxj);
        CBVtrue        = mean(CBVtrue(:));
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
    REMap = abs(CBVest-CBV)./CBV;

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
    title('Relative Error in CBV |CBVest-CBVtrue|/CBVtrue')

end



%% write image to results

if writeData
    
    
    
    basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
    
    path  = './graphics/';
    fbase = ['createindicator',indicatorcalc,'-',basenameindicator];
    fext  = ['-spacered-' num2str(blockSize(1)) '-perf.mat'];
    fname = [path,fbase,fext];
    
    %do the job
    if ~exist(fname,'file')
        save(fname,'CBV','CBVest');
    else
        save(fname,'CBV','CBVest','-append');
    end
    

    
end



%%

if showPartition
    figure(4);clf;
    set(4,'name','Partition of the image');
    imagesc(partitionImage);
    axis image;
end