% ------------------------------------------------------------------------- 
%                    E12f_SimulateConcatinatedFlow
% 
% For a voxel in the top-row, a residue function 
%
%   I  := I1\conv...\conv Il is set up for
%   Ii := Fi*exp(-Fi*t) and Fi local flow constant
%
% Then Ci is setup by convolving I with aif and compared to Cmat(i).
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
l             = 10;

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

I = zeros(kH,l);
for i = 1:l

    %setup voxel flow
    q1    = qmat{1}(1,i);
    q2    = qmat{2}(1,i);
    F     = (q1 + q2); %flow in mm^3/s
    FTrue = (F + abs(Fmat(1,i)))/hd;

    %setup voxel porosity
    phi = phimat(1,i);
    
    %setup Ii
    G      = FTrue/phi;
    I(:,i) = G*exp(-G*timelineH(:));
    
end

%% simulate C

%ONE IR
IR = I(:,1);
for i = 2:l
    IR = conv(IR,I(:,i))*dtH;
    IR = IR(1:kH);
end
cSim = conv(IR,aifH)*dtH;
CSim = phi*cSim(1:kH);

%subsequent convolution
ci = conv(IR,aifH)*dtH;
ci = ci(1:kH);
for i = 2:l
    ci = conv(IR,ci)*dtH;
    ci = ci(1:kH);
end
CSim2 = ci*phi;


%% downsample curves
nSamp = round(dt/dtH);
CSim  = CSim(1:nSamp:end);
CSim2 = CSim2(1:nSamp:end);
IR    = IR(1:nSamp:end);


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
    
    subplot(1,3,1);
    plot(timeline,CData,timeline,CRec);
    legend('CDisc','CRec');    
    title('C');
    
    subplot(1,3,2);
    plot(timeline,IRec,timeline,IR);
    legend('IRec','IR');
    ti = sprintf('FRec=%1.4f, FTrue=%1.4f',FRec,FTrue);
    title(ti);
    
    subplot(1,3,3);
    plot(timeline,IRec/norm(IRec),timeline,IR/norm(IR));
    legend('IRec/|IRec|','IR/|IR|');

end


if showC
    figure(1);clf;
    subplot(1,2,1);
    plot(timeline,CData,timeline,CSim,timeline,CSim2);
    legend('CData','C-IR','C-Sub');
    subplot(1,2,2);
    plot(timeline,CData,timeline,CSim2);
    legend('CData','C-Sub');    
    
end

