% ------------------------------------------------------------------------- 
%                               E12d_DeconvolutionOnLargerRegion 
% 
% Loads PDE model, averages a region around one voxel and does a 
% deconvolution on the mean uptake curve.
%
% CData = mean(Ci) where Ci is the average of concentrations in the grid
% around voxel idxij with given width wi and height he
%
% |__|__|__|__|__|
% |__|__|__|__|__|
% |__|__|__|__|__|   
% |__|__|ij|__|__|
% |__|__|__|__|__|
% |__|__|__|__|__|
% |__|__|__|__|__|
%
%  In the above example wi = 2 and he = 3
%
% 
%                                      (c)Constantin Sandmann, 27-Feb-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clear;
clc;
close all;

load smallDataSet

%% main variables

idx  = [32,32]; %center voxel
wi   = 3;
he   = 3;
OI   = 2e-4;

%% some basic variables
mk = size(Cmat);
m  = mk(1:3);
k  = mk(end);
dt = timeline(2)-timeline(1);

%setup indices
idxi = (idx(1)-he:idx(1)+he);
idxj = (idx(2)-wi:idx(2)+wi);

topRow = idx(1)-he;
botRow = idx(1)+he; 

rightCol = idx(2)+wi;
leftCol  = idx(2)-wi;


%% average input curve

%setup mask
mask              = false(m);
mask(idxi,idxj,1) = true;

%setup CData
CData = reshape(Cmat,[],k);
CData = mean(CData(mask(:),:),1)';


%% get flow

%get inflow from top
FTop   = qmat{1}(topRow,idxj,1);
%get inflow from left
FLeft   = qmat{1}(idxi,leftCol,1);
%just to be sure: Get outflow at bottom
FBot   = qmat{1}(botRow,idxj,1);
%just to be sure: Get outflow at right
FRight   = qmat{1}(idxi,rightCol,1);

%total flow
FIn  = sum(FTop) + sum(FLeft);
FOut = sum(FBot) + sum(FRight);

%perfusion
hd   = prod(prm.h);
vol  = nnz(mask)*hd;
FRef = FIn/vol;


fprintf('FIn=%1.4f, FOut=%1.4f, (F-FOut)/F = %1.4f\n',FIn,FOut,(FIn-FOut)/FIn);

%% get AIF

%get concentrations
c    = bsxfun(@rdivide,Cmat,phimat);
cTop  = squeeze(c(topRow,idxj,1,:));
cLeft = squeeze(c(idxi,leftCol,1,:));

%get aifs
aifTop  = sum(bsxfun(@times,cTop,FTop(:)),1);
aifLeft = sum(bsxfun(@times,cLeft,FLeft(:)),1);
aif     = (aifTop + aifLeft)./FIn;


%% deconvolution

%prepare deconvolution
A       = perfusion1c.getLinearConvolutionMatrix(aif,dt);
[U,S,V] = svd(A);

% do deconvolution
[FRec,IRec,CRec] = perfusion1c.linearDeconvolution(CData,timeline,OI,U,S,V);

%% show results
figure(2);clf;
set(2,'name','Deconvolution Results');

subplot(1,2,1);
plot(timeline,CData,timeline,CRec);
legend('CData','CRec');    
title('C');

subplot(1,2,2);
plot(timeline,IRec);
legend('IRec');
ti = sprintf('FRec=%1.4f, FRef=%1.4f',FRec,FRef);
title(ti);
