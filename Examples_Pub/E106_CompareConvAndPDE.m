% ------------------------------------------------------------------------- 
%                      E108_CompareConvAndPDE
% 
% For perfusion parameters taken from a specific voxel of the PDE-Model,
% a 1D uptake curve is simulated using
% 
% 2) Analytic: Convolution
%       C(t) = F_Ref*\int_0^t aif(t)*exp(-F_Ref/phi*(t-s)) ds 
% where F_Ref = F/hd;
%
% This curve is then compared to the actual curve from the PDE-Model.
%
% 
% 
%                                      (c)Constantin Sandmann, 09-Mar-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;
close all;

%% setup variables from software phantom

load smallDataSet

%setup timeline
timelineH = linspace(0,90,1e5);

%setup other variables
kH   = numel(timelineH);
hd   = prod(prm.h);   %voxel volume in mm^3
dtH  = timelineH(2)-timelineH(1);
dt   = timeline(2)-timeline(1);

%setup idx of voxel where the deconvolution is to be performed
idx = [32,35];

%setup voxel flow
q1  = qmat{1}(idx(1),idx(2));
q2  = qmat{2}(idx(1),idx(2));
F   = (q1 + q2); %flow in mm^3/s

%setup voxel CBV
phi = phimat(idx(1),idx(2));

%setup AIF
c1   = squeeze(Cmat(idx(1),idx(2)-1,:))./phimat(idx(1),idx(2)-1);
c2   = squeeze(Cmat(idx(1)-1,idx(2),:))./phimat(idx(1)-1,idx(2));
aif  = (q1.*c1 + q2.*c2)./(q1+q2); %weighting by amount of flow (q1 and q2)
aif  = interp1(timeline,aif,timelineH);

%% analytic solution from data
FTrue = F/hd;

%setup impuls-response function I
I = FTrue*exp(-FTrue/phi*timelineH);

%setup convolution
C    = conv(I,aif)*dtH;
CAna = C(1:kH)';

%% solution from data
%get single tissue curve
CData = squeeze(Cmat(idx(1),idx(2),:));

%% downsampling to data-grid
step        = round(dt/dtH);
timelineMod = timelineH(1:step:end);
CAna        = CAna(1:step:end);


%% show results

figure(1);clf;
set(1,'name','Discrete vs. Analytic');
plot(timeline,CAna,timeline,CData);
legend('C_{Ana}','C_{data}');

return;
%% save results for paper

figure(1);clf;
plot(timeline,CAna,timeline,CData,'linewidth',3);
set(gca,'fontsize',20)
legend('Convolution','PDE');
ylabel('concentration [mmol/mm^3]');
xlabel('time [s]');
export_fig ConvVsPDE.pdf -transparent