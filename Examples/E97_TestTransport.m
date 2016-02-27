% ------------------------------------------------------------------------- 
%                               E98_TestPressure
% 
% Checks, if p fulfills the relationship -K/mu \Delta p = Q.
% To do this, the laplacian is calculated
% 
% 
% 
%                                          (c)Constantin Heck, 06-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 



%% prepare data
%{
clear
close all;
clc


% load data
E = load('largeDataset.mat');

%get physiological dimensions
prm      = E.prm;
m        = prm.dim(1:2);
h        = prm.h;
hd       = prod(h);
voxelvol = prod(prm.h);
dim      = numel(h);
timeline = E.timeline;

%get C, q and phi
C    = E.Cmat;
AIF  = reshape(E.aifval,[],1);
phi  = E.phimat;
F    = E.Fmat;
q    = E.qmat;
k    = numel(AIF);
res  = zeros([m,k]);
%}

%% initialize variables
dt       = 0.001;
timeline = (0:dt:90);
k        = numel(timeline);
AIF      = perfusion1c.createaif('gamma',timeline);


%% do the job
% check if dC/dt = \int_{\partial Omega_i} q*c dx - Fmat*c + Fmat*aif

c = zeros([m,k]);
C = zeros([m,k]);


for i = 1:k-1
    
    %initialize variables: C and c
    Ci     = squeeze(C(:,:,i));
    ci     = squeeze(c(:,:,i));
    AIFi   = AIF(i);
    deltaT = timeline(i+1)-timeline(i);
    
    %fluxes
    fluxin  = zeros(m);
    fluxout = zeros(m);
    
    %INFLOW Top->Bottom
    cc     = [zeros(1,m(2));ci];
    cin    = cc.*q{1};
    fluxin = fluxin + cin(1:end-1,:);
    
    %INFLOW Left->Right
    cc     = [zeros(m(1),1),ci];
    cin    = cc.*q{2};
    fluxin = fluxin + cin(:,1:end-1);
    
    %OUTFLOW Top->Bottom
    cc      = [zeros(1,m(2));ci];
    cout    = cc.*q{1};
    fluxout = fluxout + cout(2:end,:);
    
    %OUTFLOW Left->Right
    cc      = [zeros(m(1),1),ci];
    cout    = cc.*q{2};
    fluxout = fluxout + cout(:,2:end);
       
    %sources
    sources = F;
    sources(sources<0) = 0;
    sources = sources*AIFi;
    
    %sinks
    sinks = F;
    sinks(sinks>0) = 0;
    sinks = -sinks.*ci;
    
    %setup update
    C(:,:,i+1) = Ci + deltaT*(fluxin - fluxout + sources - sinks)./hd;
    c(:,:,i+1) = C(:,:,i+1)./phi;
    
    
end

%% plot results
tSamp = 1;
step  = round(tSamp/deltaT);

Csmall = C(:,:,1:step:k);
ksmall = size(Csmall,3);
ca     = [min(Csmall(:)),max(Csmall(:))];

figure(1);clf;
for i = 1:ksmall
    Ci = Csmall(:,:,i);
    imagesc(Ci);
    colormap gray;
    axis image;
    caxis(ca);
    pause(.5);
end


return;

%% get errors
% tSamp  = 0.001;
% step   = round(tSamp/deltaT);
% C      = C(:,:,1:step:k);

CRef = squeeze(E.Cmat);
RE   = norm(C(:)-CRef(:))/norm(CRef(:));
RE