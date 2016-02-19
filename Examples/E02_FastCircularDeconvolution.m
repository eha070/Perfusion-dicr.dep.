% ------------------------------------------------------------------------- 
%                    E02_CircularDeconvolutionConvolution
% 
% Circular deconvolution on one of the synthetic flow models.
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



%setup the flow to compare the results to
trueFlow = 'perfusion';
% trueFlow = 'averaging'; %do NOT use this for indicatorcalc='conv'



%which results to show?
showFlowMaps       = 1;
showMultipleCurves = 0; %remember to setup indices idxiD and idxjD
showSingleCurve    = 1;
writeImage         = 0;
saveData           = 0;
saveSingleCurve    = 0;



%setup threshold
thres = 1e-1; %probably needs to be tuned



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
timelineHigh = linspace(timelineLow(1),2*timelineLow(end),2*k);
AIFlow       = AIF(idxsamp);

klow    = numel(idxsamp);
mklow   = [m,klow];


%% deconvolution

[CBFest,phi,Iest,Cest,para] = perfusion1c.fastPerfusionAnalysis(Clow,AIFlow,timelineLow,thres);

%reshape
CBFest = reshape(CBFest,m);
Iest   = reshape(Iest,[m,2*klow]);
Cest   = reshape(Cest,[m,2*klow]);




%% show results


if showFlowMaps
    %setup "segmented" Flow CBF

    %setup scaling
    tmp  = 100*60*[reshape(CBF,[],1);reshape(CBFest,[],1)];
    cmin = min(tmp);
    cmax = max(tmp);

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

end


%% show some curves


if showMultipleCurves

    %setup curves to display
    idxiD = (1:30:64);
    idxjD = (1:30:64);    
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


    REF = (CBFest(idxiD,idxjD)-CBF(idxiD,idxjD))./CBF(idxiD,idxjD);
    REF = mean(REF(:));


    figure(2);clf;

    subplot(1,2,1);
    plot(timelineHigh,Irec,'-r','lineWidth',3);
    title(sprintf('Reconstructed I, Average RE in Flow=%1.4f',REF))

    subplot(1,2,2);
    plot(timelineLow,Ctr,'-b',timelineHigh,Crec,'-r','lineWidth',3);
    title(sprintf('True C'))
    
end

%% save single curve


if saveSingleCurve
    
    
    %setup squeeze function
    s = @(v) squeeze(v(:));    
    
    %position of single curve
    pos = [50,50];
%     pos = randi(m(1),2,1);
    
    %get the stuff
    corig = s(Clow(pos(1),pos(2),:));
    cest  = s(Cest(pos(1),pos(2),:));
    iest  = s(s(Iest(pos(1),pos(2),:)));
    
    %double-check
    tau   = timelineLow(2)-timelineLow(1);
    ctest = conv(AIFlow,iest);
    ctest = tau*ctest(1:k);


    figure(3);clf;
    plot(timelineLow,corig,timelineLow,ctest,'LineWidth',3);
    legend('Measured C','Model Approximation of C');
    xlabel('Time (s)')
    ylabel('Concentration (mol/l)')
    set(gca,'FontSize',15)

    export_fig ./figs/C-and-Crec-PDE.eps -transparent

    pause(1);
    
    figure(3);clf;
    plot(timelineLow,iest(1:k),'LineWidth',3);
    xlabel('Time (s)')
    legend('I')
    set(gca,'FontSize',15)

    export_fig ./figs/Irec-PDE.eps -transparent    
    
end





%% results for a single voxel

if showSingleCurve
    
    %position of single curve
    pos = [32,32];

    s = @(v) squeeze(v(:));

    figure(3);clf;
    set(3,'name','Comparison on single voxel');

    subplot(1,2,1);
    plot(timelineLow,s(Clow(pos(1),pos(2),:)),timelineHigh,s(Cest(pos(1),pos(2),:)),'LineWidth',3);
    title('Ctrue (blue) and reconstructed C (red) at position pos')
    
    subplot(1,2,2);
    plot(timelineHigh,s(Iest(pos(1),pos(2),:)),'LineWidth',3);
    title('estimated I at position pos')
end


%% testing
%{
C = s(Clow(pos(1),pos(2),:));
[F,~,Iest,Cest,para] = perfusion1c.fastPerfusionAnalysis(C,AIFlow,timelineLow,thres,'mode','bcSVD');
plot(timelineLow,C,timelineHigh,Cest)
% plot(timelineLow,AIFlow)
% plot(timelineLow,Iest)
%}



%% results for a single voxel

if writeImage

    %setup scaling
    tmp  = 100*60*[CBF(:);reshape(CBFest(CBFest~=0),[],1)];

    %setup Relative Error
    REMap = 100*abs(CBFest-CBF)./CBF;
    
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
    fpicname = ['./figs/','recCirc-',fbase,'-RE-Flow.eps'];
    eval(['export_fig ',fpicname,' -transparent']);
   
    
end



%% save the data

if saveData
    
    %setup base filename
    basenameindicator = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T');
    fbase = ['createindicator',indicatorcalc,'-',basenameindicator];
    
    %save files
    fmatname = ['./results/','recCirc-',fbase,'.mat'];
    
    fprintf('Saving Results...');
    tic;
    save(fmatname);
    fprintf('...done. Elapsed time: %1.2fs \n',toc);
    
end
