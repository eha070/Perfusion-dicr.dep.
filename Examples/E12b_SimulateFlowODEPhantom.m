% ------------------------------------------------------------------------- 
%                               E12b_SimulateFlowODEPhantom
% 
% For perfusion parameters taken from a specific voxel of the PDE-Model,
% a 1D uptake curve is simulated by using
% 
% 1) Discrete: Forward-Euler
%       C'(t) = F*(aif(t) - c(t))
%
% 2) Analytic: Convolution
%       C(t) = F_Ref*\int_0^t aif(t)*exp(-F_Ref/phi*(t-s)) ds 
% where F_Ref = F/hd;
%
% These two curves are then compared to the actual curve from the
% PDE-Model.
%
% 
% 
%                                      (c)Constantin Sandmann, 27-Feb-2016 
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
idx = [32,32];

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

%% solution by discretization

c = zeros(kH,1);
C = zeros(kH,1);
for i = 1:kH-1
    
    incr = F/hd*(aif(i) - c(i));
    
    C(i+1) = C(i) + dtH*incr;
    c(i+1) = C(i+1)/phi;
    
end

CDisc = C;
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
CDisc       = CDisc(1:step:end);


%% show results

figure(1);clf;
set(1,'name','Discrete vs. Analytic');
plot(timeline,CAna,timeline,CDisc,timeline,CData);
legend('C_{Ana}','C_{Disc}','C_{data}');
