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
% indicatorcalc = 'conv';
indicatorcalc = 'PDE';



%setup the flow to compare the results to
trueFlow = 'perfusion';
% trueFlow = 'averaging'; %do NOT use this for indicatorcalc='conv'



%which results to show?
showFlowMaps       = 1;
showMultipleCurves = 0; %remember to setup indices idxiD and idxjD
showSingleCurve    = 0;
writeData          = 0;
showPartition      = 0;


%setup size of blocks to run the deconvolution on
blockSize = 1*[1,1];


%setup oscillation index OI
% OI = .001;
OI = .01;


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

switch trueFlow
    case 'perfusion'
        CBF = E.perfmat; 
        CBF(1,1) = 2*CBF(1,1);
        CBF(end,end) = 2*CBF(end,end);
        
    case 'averaging'
        qcc = perfusion1c.convertFlowStagToCC(E.qmat);
        CBF = .5*(qcc{1} + qcc{2});

    otherwise
        error('Value of trueFlow not right');
end


%clear memory
clearvars E D;


% return;

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


%% setup area where to run the deconvolution


%initialize variables to store results
CBFest = zeros(m);
Iest   = zeros([m,2*klow]);
Cest   = zeros([m,2*klow]);
Cblock = zeros([m,klow]);

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
        idxi(idxi>csi(1)) = [];
        idxj(idxj>csi(2)) = [];
        
        %get current blockSize
        cBlockSize = [numel(idxi),numel(idxj)];
        cBlockn    = prod(cBlockSize);
       
        
        %create image to display the blocks
        partitionImage(idxi,idxj) = (i-1)*nBlocksj + j;

        
        %do the deconvolution
        Cij = Clow(idxi,idxj,:);
        Cij = reshape(Cij,[],k);
        Cij = mean(Cij,1)';
        [F,Irec,Crec] = perfusion1c.circularDeconvolution(Cij,timelinelow,OI,U,S,V);
        
        
        %prepare Irec and Crec to store them
        Irec = reshape(Irec,1,1,[]);
        Irec = repmat(Irec,cBlockSize(1),cBlockSize(2),1);
        
        Crec = reshape(Crec,1,1,[]);
        Crec = repmat(Crec,cBlockSize(1),cBlockSize(2),1);
        
        %prepare Cbloc
        Cij = reshape(Cij,1,1,[]);
        Cij = repmat(Cij,cBlockSize(1),cBlockSize(2),1);

        

        %store results
        CBFest(idxi,idxj)   = F;
        Iest(idxi,idxj,:)   = Irec;
        Cest(idxi,idxj,:)   = Crec;
        Cblock(idxi,idxj,:) = Cij;

        %fill the waitbar
        perc = (nBlocksj*(i-1) + j)/nBlocks;
        waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));       
        
        
        %prepare true flow to block-size
        Ftrue          = CBF(idxi,idxj);
        Ftrue          = mean(Ftrue(:));
        CBF(idxi,idxj) = Ftrue;
       
       
   end
end

delete(h);


%% showFlowMaps

if showFlowMaps
    
    %setup scaling cmax
    cScale = 1;
    
    %setup color axis
    tmp  = [CBF(:);CBFest(:)]*100*60;
    cmin = min(tmp);
    cmax = cScale*max(tmp);

    %setup 
    REMap = abs(CBFest-CBF)./CBF;

    %show flow
    figure(1);clf;
    colormap jet(512);
    set(gcf,'name','Results for Convolution Model')


    subplot(1,3,1);
    imagesc(CBF*100*60);
    caxis([cmin,cmax]);
    axis image;
    title('True flow (ml/min/100ml)')

    subplot(1,3,2);
    imagesc(CBFest*100*60);
    caxis([cmin,cmax]);
    axis image;
    title('Estimated flow (ml/min/100ml)')

    %show division
    subplot(1,3,3);
    imagesc(REMap);
    caxis([0,1]);
    axis image;
    title('Relative Error in Flows |Fest-Ftrue|/Ftrue')
    
    medianRE = median(REMap(:));
    fprintf('medianRE: %1.4f\n',medianRE);

end

%% show some curves

if showMultipleCurves

    %setup curves to display
    idxiD = (1:10:64);
    idxjD = (1:10:64);    
    numi   = numel(idxiD)*numel(idxjD);


    %reshape everything
    Ctr = Cblock(idxiD,idxjD,:);
    Ctr = reshape(Ctr,numi,[]);
    Ctr = Ctr';

    Crec = Cest(idxiD,idxjD,:);
    Crec = reshape(Crec,numi,[]);
    Crec = Crec';

    Irec = Iest(idxiD,idxjD,:);
    Irec = reshape(Irec,numi,[]);
    Irec = Irec';


    REF = (CBFest(idxiD,idxjD)-CBF(idxiD,idxjD))./CBF(idxiD,idxjD);
    REF = mean(REF(:));


    figure(2);clf;

    subplot(1,2,1);
    plot(timelinelow2,Irec);
    title(sprintf('Reconstructed I, Average RE in Flow=%1.4f',REF))

    subplot(1,2,2);
    plot(timelinelow,Ctr,'-b',timelinelow2,Crec,'-r');
    title(sprintf('True C in blue, reconstructed C in red'))

end

%% results for a single voxel

if showSingleCurve
    
    %setup squeeze function
    s = @(v) squeeze(v(:));    
    
    %position of single curve
    pos = [50,50];
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





%% write image to results

if writeData
    
    
    
    basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
    
    path  = './graphics/';
    fbase = ['createindicator',indicatorcalc,'-',basenameindicator];
    fext  = ['-spacered-' num2str(blockSize(1)) '-perf.mat'];
    fname = [path,fbase,fext];
    
    %do the job
    CBFCirc = CBFest;
    if ~exist(fname,'file')
        save(fname,'CBF','CBFCirc');
    else
        save(fname,'CBFCirc','CBF','-append');
    end
    

    
end


%%

saveSingleCurve = 0;
if saveSingleCurve
    
    
    %setup squeeze function
    s = @(v) squeeze(v(:));    
    
    %position of single curve
    pos = [50,50];
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
    plot(timelinelow,corig,timelinelow,ctest,'LineWidth',3);
    legend('Measured C','Model Approximation of C');
    xlabel('Time (s)')
    ylabel('Concentration (mol/l)')
    set(gca,'FontSize',15)

    export_fig ./figs/C-and-Crec-PDE.eps -transparent
    
    figure(3);clf;
    plot(timelinelow,iest(1:k),'LineWidth',3);
    xlabel('Time (s)')
%     ylabel('Concentration (mmol/mm^3)')
    legend('I')
    set(gca,'FontSize',15)

    export_fig ./figs/Irec-PDE.eps -transparent    
    
end




%%

if showPartition
    figure(4);clf;
    set(4,'name','Partition of the image');
    imagesc(partitionImage);
    axis image;
end