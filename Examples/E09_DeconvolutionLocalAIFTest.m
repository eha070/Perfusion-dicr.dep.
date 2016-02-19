% ------------------------------------------------------------------------- 
%                    E09_DeconvolutionLocalAIFTest
% 
% Do a deconvolution for one voxel with upstream voxels as AIF.
% The following model is assumed:
%
%  Ci = CBFi (R * AIF)
%
% where * denotes convolution, CBFi and R are to be determined. 
% The AIF is chosen as 
%
%   AIF = (C1/phi1+C2/phi2)/2
%
% where phi denotes the porosity and 
%
%       |  C1   |
%    ---|---|---|----
%       |   v   |
%    C2 -> Ci   |
%       |       |
%    ---|-------|----
%       |       |
% 
% The ground-truth value should hence be given by
%
%  CBFi = (q1 + q2)/(2*hd)
%
% here qj denotes the (absolute) flow from Cj to Ci.
% Note that a division by hd is performed to convert the (absolute) flow 
% [mm^3/s] to perfusion [mm^3/s/mm^3].
% 
%
%  
%                                          (c)Constantin Heck, 19-Nov-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clc;
close all;



%setup oscillation index OI
OI = .0002; %for PDE;
% OI = .1;

%setup idx of voxel where the deconvolution is to be performed
idx = [35,35];



%% prepare data

%get single tissue curve
CLoc    = squeeze(Cmat(idx(1),idx(2),:));

%get aif
aifLoci = squeeze(Cmat(idx(1),idx(2)-1,:))./phimat(idx(1),idx(2));
aifLocj = squeeze(Cmat(idx(1)-1,idx(2),:))./phimat(idx(1),idx(2));
aifLoc  = (aifLoci + aifLocj)./2;


%% setup area where to run the deconvolution

%get deconvolution matrix
deltaT  = timeline(2)-timeline(1);
A       = perfusion1c.getLinearConvolutionMatrix(aifLoc,deltaT);
fprintf('Starting SVD...');
tic; [U,S,V] = svd(A);
fprintf('...done. Elapsed time: %1.3fs\n',toc);


%do the deconvolution
[CBFrec,Irec,Crec] = perfusion1c.linearDeconvolution(CLoc,timeline,OI,U,S,V);


%% show results

%show flow
figure(1);clf;

subplot(1,3,1);
plot(timeline,aifLoc,timeline,CLoc,'lineWidth',3);
legend('aifLoc','CLoc')
title(ti);

subplot(1,3,2);
plot(timeline,Irec,'lineWidth',3);
ti = sprintf('Impuls-Response Function: CBF=%1.4f',CBFrec*100*60);
title(ti);

subplot(1,3,3);
plot(timeline,CLoc,timeline,Crec,'lineWidth',3);
legend('true C','recovered C');
ti = sprintf('Reconstruction');
title(ti);



%% setup true flow

%qmat is given absolute with units [mm^3/s].
%It needs to be converted to perfusion [mm^3/s/mm^3]
%divide by voxel-volume
CBFtr = (qmat{1}(idx(1),idx(2)) + qmat{2}(idx(1),idx(2)))./(2*hd);
CBVtr = CBV(idx(1),idx(2));


%% median errors

RECirc = abs(CBFrec-CBFtr)./CBFtr*100;
REMS   = abs(CBFrecMS-CBFtr)./CBFtr*100;
RECBV  = abs(CBVrec-CBVtr)./CBVtr*100;

fprintf('RE in Lin: \t RE=%1.2f%% \n',RECirc);
fprintf('RE in MS: \t RE=%1.2f%% \n',REMS);
fprintf('RE in CBV: \t RE=%1.2e%% \n',RECBV);


%% 
saveImage = 0;

if saveImage
    
    figure(2);clf;
    plot(timeline,Irec,'lineWidth',3);
    xlabel('Time (s)')
%     ylabel('Concentration (mmol/mm^3)')
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
    
    
    
