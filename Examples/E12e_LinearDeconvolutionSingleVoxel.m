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
aifMode       = 'global';
idx           = [31,1];

OI            = 1e-2;

%% setup main variables

load smallDataSet

%setup timeline
k        = numel(timeline);
dt       = timeline(2)-timeline(1);
hd       = prod(prm.h);   %voxel volume in mm^3

%% Ground-Truth Flow

%setup voxel flow
q1    = qmat{1}(idx(1),idx(2));
q2    = qmat{2}(idx(1),idx(2));
F     = (q1 + q2); %flow in mm^3/s
FTrue = (F + abs(Fmat(idx(1),idx(2))))/hd;

%setup voxel porosity
phi = phimat(idx(1),idx(2));



%% get data curve
CData = squeeze(Cmat(idx(1),idx(2),1,:));


%% get AIF for deconvolution

switch aifMode
    case 'local'
        c1  = squeeze(Cmat(idx(1),idx(2)-1,:))./phimat(idx(1),idx(2)-1);
        c2  = squeeze(Cmat(idx(1)-1,idx(2),:))./phimat(idx(1)-1,idx(2));
        aif = (q1.*c1 + q2.*c2)./(q1+q2); %weighting by amount of flow (q1 and q2)
    case 'global'
        aif  = aifval;
end


%% recover flow by deconvolution

%prepare deconvolution
A       = perfusion1c.getLinearConvolutionMatrix(aif,dt);
tic; fprintf('Starting SVD...');
[U,S,V] = svd(A);
fprintf('finished: %1.2fs\n',toc);

%do deconvolution
[FRec,IRec,CRec] = perfusion1c.linearDeconvolution(CData,timeline,OI,U,S,V);

%% show results

if showDeconv
    figure(2);clf;
    set(2,'name','Deconvolution Results');
    
    subplot(1,2,1);
    plot(timeline,CData,timeline,CRec);
    legend('CDisc','CRec');    
    title('C');
    
    subplot(1,2,2);
    plot(timeline,IRec);
    legend('IRec');
    ti = sprintf('FRec=%1.4f, FTrue=%1.4f',FRec,FTrue);
    title(ti);

end