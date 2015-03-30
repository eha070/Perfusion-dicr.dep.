% ------------------------------------------------------------------------- 
%                    E04_OneBigCompartment
% 
% Circular deconvolution on the synthetic flow model, set up with
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
indicatorcalc = 'conv';


%setup oscillation index OI
OI = .0005; %for PDE;
% OI = .001; %for CONV;


%prepare downsampling of data (shorter runtime for SVD)
% It will hold: Clow = Clow(:,:,1,1:step:end);
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


%average to one single uptake curve Cav
Cav = reshape(Clow,n,k);
Cav = mean(Cav,1)';







%% setup area where to run the deconvolution

%get deconvolution matrix
deltaT  = timelinelow(2)-timelinelow(1);
A       = perfusion1c.getLinearConvolutionMatrix(AIFlow,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);


%do the deconvolution
[CBFrec,Irec,Crec] = perfusion1c.linearDeconvolution(Cav,timelinelow,OI,U,S,V);

%do the maximum slope technique
CBFrecMS = perfusion1c.maximumSlope(Cav,timelinelow,AIFlow);

%get CBV
CBVrec = perfusion1c.cbvEstimation(Cav,timelinelow,AIFlow);


%% show results

%show flow
figure(1);clf;
subplot(1,2,1);
plot(timelinelow,Cav,timelinelow,Crec,'lineWidth',3);
legend('true C','recovered C');
ti = sprintf('Maximum-Slope: CBF=%1.4f',CBFrecMS*100*60);
title(ti);

subplot(1,2,2);
plot(timelinelow,Irec,'lineWidth',3);
ti = sprintf('Impuls-Response Function: CBF=%1.4f',CBFrec*100*60);
title(ti);





%% median errors
CBFtr = mean(CBF(:));
CBVtr = mean(CBV(:));

RECirc = abs(CBFrec-CBFtr)./CBFtr*100;
REMS   = abs(CBFrecMS-CBFtr)./CBFtr*100;
RECBV  = abs(CBVrec-CBVtr)./CBVtr*100;

fprintf('RE in Circ: \t RE=%1.2f%% \n',RECirc);
fprintf('RE in MS: \t RE=%1.2f%% \n',REMS);
fprintf('RE in CBV: \t RE=%1.2e%% \n',RECBV);


%% 
saveImage = 1;

if saveImage
    
    figure(2);clf;
    plot(timelinelow,Irec,'lineWidth',3);
    xlabel('Time (s)')
%     ylabel('Concentration (mmol/mm^3)')
    legend('I')
    set(gca,'FontSize',15)
    
    export_fig ./figs/Irec-conv.eps -transparent
    

    figure(3);clf;
    plot(timelinelow,Cav,timelinelow,Crec,'lineWidth',3);
    xlabel('Time (s)')
    ylabel('Concentration (mmol/mm^3)')
    legend('C','Model Approximation of C')
    set(gca,'FontSize',15)

    export_fig ./figs/C-and-Crec-conv.eps -transparent
     
end
    
    
    
