% ------------------------------------------------------------------------- 
%                    E03_LinearDeconvolution
% 
% Cirvular deconvolution on the synthetic flow model, set up with
% exponential residue functions
% 
% 
% 
%                                          (c)Constantin Heck, 19-Nov-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clear;
clc;
close all;


%setup which flow-calculation to use
% indicatorcalc = 'conv';
indicatorcalc = 'PDE';


%setup area where to run the deconvolution
idxi  = (1:64);
idxj  = (1:64);


%setup oscillation index OI
OI = .5; %determined experimentally;


%prepare downsampling of data (shorter runtime for SVD)
% It will hold: Clow = Clow(:,:,1,1:step:k);
step = 1;




%% load data

%settings
prm               = settings;
basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
basenameFlow      = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');



%setup paths
pathloadFlow = ['synt-createflowTPFA-' basenameFlow '.mat'];
switch indicatorcalc
    case 'PDE'
        pathload = ['synt-createindicatorpde-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
    case 'conv'
        pathload = ['syntconv-createindicatorconv-' basenameindicator '-red-' int2str(prm.stepred) '.mat'];
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
CBF = E.perfmat;

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


%% setup area where to run the deconvolution

%downsmaple the data
numi        = numel(idxi)*numel(idxj);


%initialize variables to store results
CBFest = zeros(m);
Iest   = zeros([m,klow]);
Cest   = zeros([m,klow]);

%get deconvolution matrix
deltaT  = timelinelow(2)-timelinelow(1);
A       = perfusion1c.getLinearConvolutionMatrix(AIFlow,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);

%do the job voxelwise
h = waitbar(0);
for i = idxi;
    for j = idxj;
                    
       
        %do the deconvolution
        Cij = squeeze(Clow(i,j,:));
        [F,Irec,Crec] = perfusion1c.linearDeconvolution(Cij,timelinelow,OI,U,S,V);
        
        %store results
        CBFest(i,j) = F;
        Iest(i,j,:) = Irec;
        Cest(i,j,:) = Crec;
        
        %fill the waitbar
        perc = (numel(idxj)*(i-min(idxi)) + j-min(idxj)+1 )/numi;
        waitbar(perc,h,sprintf('Calculating...%2.2f%%...',perc*100));

        
    end
end

close(h);

%% show results

%setup "segmented" Flow CBFseg
seg = zeros(m);
seg(idxi,idxj) = 1;
CBFseg = CBF.*seg;

%setup scaling
tmp  = 100*60*[reshape(CBFseg(CBFseg~=0),[],1);reshape(CBFest(CBFest~=0),[],1)];
cmin = min(tmp);
cmax = max(tmp);

%setup 
REMap = abs(CBFest-CBFseg)./CBFseg;

%show flow
figure(1);clf;
colormap jet(512);
set(gcf,'name','Results for Convolution Model')


subplot(1,3,1);
imagesc(CBFseg*100*60);
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


%% show some curves

%setup indices to display using steps
stepi = 25;
stepj = 25;
idxiD  = (idxi(1):stepi:idxi(end));
idxjD  = (idxj(1):stepj:idxj(end));

%manual setup of indices to display
idxiD = (10:20:60);
idxjD = (10:20:60);

numi   = numel(idxiD)*numel(idxjD);




Ctr = Clow(idxiD,idxjD,:);
Ctr = reshape(Ctr,numi,[]);
Ctr = Ctr';

Crec = Cest(idxiD,idxjD,:);
Crec = reshape(Crec,numi,[]);
Crec = Crec';

Irec = Iest(idxiD,idxjD,:);
Irec = reshape(Irec,numi,[]);
Irec = Irec';


REF = (CBFest(idxi,idxj)-CBF(idxi,idxj))./CBF(idxi,idxj);
REF = mean(REF(:));


figure(2);clf;

subplot(1,2,1);
plot(timelinelow,Irec);
title(sprintf('Reconstructed I, Average RE in Flow=%1.4f',REF))

subplot(1,2,2);
plot(timelinelow,Ctr,'-b',timelinelow,Crec,'-r');
title(sprintf('True C'))



%%
figure(1);clf;
subplot(1,3,1);
plot(s(Clow(32,32,:)),'LineWidth',3); set(gca,'FontSize',25)

subplot(1,3,2);
plot(s(Cest(32,32,:)),'LineWidth',3); set(gca,'FontSize',25)

subplot(1,3,3);
plot(s(Iest(32,32,:)),'LineWidth',3); set(gca,'FontSize',25)

