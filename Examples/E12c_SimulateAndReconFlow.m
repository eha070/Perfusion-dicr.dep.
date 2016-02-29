% ------------------------------------------------------------------------- 
%                               E12_SimulateFlowODE 
% 
% For perfusion parameters taken from a specific voxel of the PDE-Model,
% a 1D uptake curve is simulated by using
% 
% 1) Analytic: Convolution
%       C(t) = F_Ref*\int_0^t aif(t)*exp(-F_Ref/phi*(t-s)) ds 
% where F_Ref = F/hd;
%
% We then try to recover the perfusion by a doing a deconvolution.
%
%                                      (c)Constantin Sandmann, 27-Feb-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;
close all;

showDeconv    = 1;

OI            = 1e-3;

%% setup main variables

load smallDataSet

%setup timeline
timelineH = linspace(0,90,1e5);
kH        = numel(timelineH);
dtH       = timelineH(2)-timelineH(1);

%setup idx of voxel where the deconvolution is to be performed
idx = [32,32];
hd  = prod(prm.h);   %voxel volume in mm^3

%setup voxel flow
q1  = qmat{1}(idx(1),idx(2));
q2  = qmat{2}(idx(1),idx(2));
F   = (q1 + q2); %flow in mm^3/s

%setup voxel CBV
phi = phimat(idx(1),idx(2));

%setup AIF
c1    = squeeze(Cmat(idx(1),idx(2)-1,:))./phimat(idx(1),idx(2)-1);
c2    = squeeze(Cmat(idx(1)-1,idx(2),:))./phimat(idx(1)-1,idx(2));
aifH  = (q1.*c1 + q2.*c2)./(q1+q2); %weighting by amount of flow (q1 and q2)
aifH  = interp1(timeline,aifH,timelineH);

%% analytic solution from data
FTrue = F/hd;

%setup impuls-response function I
IH = FTrue*exp(-FTrue/phi*timelineH);

%setup convolution
C     = conv(IH,aifH)*dtH;
CAnaH = C(1:kH)';


%% downsampling to data-grid
dt       = timeline(2)-timeline(1);
step     = round(dt/dtH);
timeline = timelineH(1:step:end);
CAna     = CAnaH(1:step:end);
aif      = aifH(1:step:end);
I        = IH(1:step:end);


%% recover flow by deconvolution

%prepare deconvolution
A       = perfusion1c.getLinearConvolutionMatrix(aif,dt);
tic; fprintf('Starting SVD...');
[U,S,V] = svd(A);
fprintf('finished: %1.2fs\n',toc);

%do deconvolution
[FRec,IRec,CRec] = perfusion1c.linearDeconvolution(CAna,timeline,OI,U,S,V);

%% show results
if showDeconv
    figure(2);clf;
    set(2,'name','Deconvolution Results');
    
    subplot(1,2,1);
    plot(timeline,CAna,timeline,CRec);
    legend('C_Disc','CRec');    
    title('C');
    
    subplot(1,2,2);
    plot(timeline,I,timeline,IRec);
    legend('I','IRec');
    ti = sprintf('FRec=%1.4f, FTrue=%1.4f',FRec,FTrue);
    title(ti);

end

