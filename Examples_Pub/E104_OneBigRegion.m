% ------------------------------------------------------------------------- 
%                    E104_OneBigRegion
% 
% Averages the complete software phantom to estimate CBF using different
% techniques.
% 
% 
% 
% 
%                                      (c)Constantin Sandmann, 24-Feb-2017 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clc;
close all;

showCurves = 1;
saveCurvesForPubliction = 1;


%setup oscillation index OI
OI = .00005;



%% load and prepare data

load smallDataSet.mat

Cav = reshape(Cmat,n,k);
Cav = mean(Cav,1)';


%% setup area where to run the deconvolution

%get deconvolution matrix
deltaT  = timeline(2)-timeline(1);
A       = perfusion1c.getLinearConvolutionMatrix(aifval,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);


%do the deconvolution
[CBFrec,Irec,Crec] = perfusion1c.linearDeconvolution(Cav,timeline,OI,U,S,V);

%do the maximum slope technique
CBFrecMS = perfusion1c.maximumSlope(Cav,timeline,aifval);

%get CBV
CBVrec = perfusion1c.cbvEstimation(Cav,timeline,aifval);

%some text output
%CBFtr = sum(Fmat(Fmat>0))/Hd;
CBVtr = 0.05;

RECirc = abs(CBFrec-CBFtr)./CBFtr*100;
REMS   = abs(CBFrecMS-CBFtr)./CBFtr*100;
RECBV  = abs(CBVrec-CBVtr)./CBVtr*100;

fprintf('RE in Circ: \t RE=%1.2f%% \n',RECirc);
fprintf('RE in MS: \t RE=%1.2f%% \n',REMS);
fprintf('RE in CBV: \t RE=%1.2e%% \n',RECBV);




%% show results
if showCurves
    %show flow
    figure(1);clf;
    subplot(1,2,1);
    plot(timeline,Cav,timeline,Crec,'lineWidth',3);
    legend('true C','recovered C');
    ti = sprintf('Maximum-Slope: CBF=%1.4f',CBFrecMS*100*60);
    title(ti);

    subplot(1,2,2);
    plot(timeline,Irec,'lineWidth',3);
    ti = sprintf('Impuls-Response Function: CBF=%1.4f',CBFrec*100*60);
    title(ti);

end

%% save images for pub

if saveCurvesForPubliction
    figure(2);clf;
    plot(timeline,Irec,'lineWidth',3);
    xlabel('Time (s)')
    legend('I')
    set(gca,'FontSize',15)

    export_fig ./figs/Irec.eps -transparent


    figure(3);clf;
    plot(timeline,Cav,timeline,Crec,'lineWidth',3);
    xlabel('Time (s)')
    ylabel('Concentration (mol/l)')
    legend('C','Model Approximation of C')
    set(gca,'FontSize',15)

    export_fig ./figs/C-and-Crec.eps -transparent

end    
    
