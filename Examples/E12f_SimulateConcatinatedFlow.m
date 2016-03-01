% ------------------------------------------------------------------------- 
%                    E12f_SimulateConcatinatedFlow
% 
% For a voxel in the top-row, a residue function 
%
%   I  := I1\conv...\conv Il is set up for
%   Ii := Fi*exp(-Fi*t) and Fi a local flow constant.
%
% Then Ci is setup by convolving I with aif and compared to Cmat(i).
% Additionally I is compared to IRec, the recovered residue function from 
% Cmat.
%
%                                      (c)Constantin Sandmann, 27-Feb-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;
close all;

showDeconv    = 1;
showC         = 1;
aifMode       = 'global';
l             = 3; %analytic solution only for these points

OI            = 1e-2;

%% setup main variables

load smallDataSet

%setup timeline
timelineH = linspace(0,90,2e6);
kH        = numel(timelineH);
dtH       = timelineH(2)-timelineH(1);

%setup other parameters
k  = numel(timeline);
dt = timeline(2)-timeline(1);
hd = prod(prm.h);   %voxel volume in mm^3

%setup AIF
aif  = aifval;
aifH = perfusion1c.getGammaAIF(timelineH/60)*1e-6;

%% setup residue functino Ii

I     = zeros(kH,l);
FTrue = zeros(l,1);
phi   = zeros(l,1);
G     = zeros(l,1);
for i = 1:l

    %setup voxel flow
    q1       = qmat{1}(1,i);
    q2       = qmat{2}(1,i);
    F        = (q1 + q2); %flow in mm^3/s
    FTrue(i) = (F + abs(Fmat(1,i)))/hd;

    %setup voxel porosity
    phi(i) = phimat(1,i);
    
    %setup Ii
    G(i)   = FTrue(i)/phi(i);
    I(:,i) = G(i)*exp(-G(i)*timelineH(:));
    
end

%% simulate C

%discrete IR
IR = I(:,1);
for i = 2:l
    IR = conv(IR,I(:,i))*dtH;
    IR = IR(1:kH);
    fprintf('%i ',i);    
end
fprintf('\n');
cSim = conv(IR,aifH)*dtH;
CSim = phi(l)*cSim(1:kH);


%analytic solution (up to now only for 2 or 3 points)
IAna2 = G(1)*G(2)/(G(1)-G(2))*(exp(-G(2)*timelineH(:)) - exp(-G(1)*timelineH(:)));
IAna  = prod(G)/((G(1)-G(2))*(G(2)-G(3))*(G(1)-G(3)))*((G(1)-G(2))*exp(-G(3)*timelineH(:)) + (G(2)-G(3))*exp(-G(1)*timelineH(:)) - (G(1)-G(3))*exp(-G(2)*timelineH(:)));

cAna = conv(IAna,aifH)*dtH;
CAna = phi(l)*cSim(1:kH);


%% downsample curves
nSamp     = round(dt/dtH);
CSim      = CSim(1:nSamp:end);
CAna      = CAna(1:nSamp:end);
IR        = IR(1:nSamp:end);
timeline2 = timelineH(1:nSamp:end);


%% get data curve
CData = squeeze(Cmat(1,l,1,:));


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
    plot(timeline,IRec,timeline2,IR,timeline2,IAna);
    legend('IRec','IR');
    ti = sprintf('FRec=%1.4f, FTrue=%1.4f',FRec,FTrue);
    title(ti);
    
end


if showC
    figure(1);clf;
    plot(timeline,CData,timeline2,CSim,timeline2,CAna);
    legend('CData','C-IR');    
end

