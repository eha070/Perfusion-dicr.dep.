% ------------------------------------------------------------------------- 
%                    E03_MaximumSlope
% 
% Maximum Slope Method on one of the synthetic flow models.
% 
% 
%                                          (c)Constantin Heck, 19-Nov-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;
close all;

%setup which flow-calculation to use
indicatorcalc = 'conv';
% indicatorcalc = 'pde';



%setup the flow to compare the results to
% trueFlow = 'perfusion';
trueFlow = 'averaging'; %do NOT use this for indicatorcalc='conv'



%which results to show?
showFlowMaps = 1;
writeImage   = 0;
saveData     = 0;


%setup area where to run the deconvolution
idxi  = (1:64);
idxj  = (1:64);


%setup oscillation index OI
OI = .01; %probably needs to be tuned



%prepare downsampling of data (shorter runtime for SVD)
%It will hold: Clow = Clow(:,:,1,1:step:k);
step = 1;



%% load data

%settings
prm               = settings;
basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
basenameFlow      = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');
foldername        = './results/';


%setup paths
pathloadFlow = [foldername,'synt-createflowTPFA-' basenameFlow '.mat'];
switch indicatorcalc
    case 'pde'
        pathload = [foldername,'synt-createindicatorpde-' basenameindicator '-red-' num2str(prm.stepred) '.mat'];
    case 'conv'
        pathload = [foldername,'syntconv-createindicatorconv-' basenameindicator '-red-' num2str(prm.stepred) '.mat'];
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
idxsamp      = (1:step:k);
Clow         = squeeze(Cmat(:,:,idxsamp));
timelineLow  = timeline(idxsamp);
timelineHigh = linspace(timeline(1),timeline(end),2*k);
AIFlow       = AIF(idxsamp);

klow    = numel(idxsamp);
mklow   = [m,klow];


%% setup area where to run the deconvolution

CBFest = perfusion1c.maximumSlope(Cmat,timelineLow,AIFlow);
CBFest = reshape(CBFest,m);
        

%% show results


if showFlowMaps
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

end








%% results for a single voxel

if writeImage

    %setup scaling
    tmp  = 100*60*[CBF(:);reshape(CBFest(CBFest~=0),[],1)];

    %setup Relative Error
    REMap = 100*abs(CBFest-CBFseg)./CBFseg;
    
    %setup figure
    figure(4);clf;
    imagesc(REMap);
    colormap jet(512);
    axis image;
    axis off;
    caxis([0,100]);
    colorbar;    
    
    
    %output RE
    medianRE = median(REMap(:));
    fprintf('Median RE: \t \t %1.2f%% \n',medianRE);
    
    %setup base filename
    basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
    fbase = ['createindicator',indicatorcalc,'-',basenameindicator];

    %save image    
    fpicname = ['./figs/','recMS-',fbase,'-RE-Flow.eps'];
    eval(['export_fig ',fpicname,' -transparent']);
        
end



%% save the data

if saveData
    
    %setup base filename
    basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
    fbase = ['createindicator',indicatorcalc,'-',basenameindicator];
    
    %save files
    fmatname = ['./results/','recMS-',fbase,'.mat'];
    
    fprintf('Saving Results...');
    tic;
    save(fmatname);
    fprintf('...done. Elapsed time: %1.2fs \n',toc);
    
end
