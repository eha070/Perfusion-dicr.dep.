% ------------------------------------------------------------------------- 
%                               E12_SimulateFlowODE 
% 
% Simulates a tissue curve C(t) = phi*c(t) based on two different ways
% 
% 1) Discrete: Forward-Euler for
%       C'(t) = F*(aif(t) - c(t))
%
% 2) Analytic: Convolution
%       C(t) = F_Ref*\int_0^t aif(t)*exp(-F_Ref/phi*(t-s)) ds 
% where F_Ref = F/hd;
%
% Then recovers F by deconvolving C with the AIF.
% 
%                                      (c)Constantin Sandmann, 27-Feb-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


showDiscVsAna = 1;
showDeconv    = 1;

%% setup main variables

%THESE PARAMETERS WORK
phi  = 0.05;  %porosity in ml/ml
F    = 50/60; %flow in ml^3/s
hd   = 100;   %voxel volume in ml
OI   = .0003;
res  = 1e-1; %sampling period
tMax = 90;
nSmp = 1e5;

% %APPROXIMATE PARAMTERS FROM EXAMPLES
% phi  = 0.05;  %porosity in ml/ml
% F    = 0.02; %flow in ml^3/s
% hd   = 0.02;   %voxel volume in ml
% OI   = 0.003;
% res  = 1e-2;
% tMax = 90;
% nSmp = 1e5;


%main input parameters
timeline = linspace(0,tMax,nSmp)';
aif      = perfusion1c.getGammaAIF(timeline/60);
dt       = timeline(2)-timeline(1);
k        = numel(timeline);


%% solution by discretization

c = zeros(k,1);
C = zeros(k,1);
for i = 1:k-1
    
    incr = F/hd*(aif(i) - c(i)); %increment in mmol/ml/s
    
    C(i+1) = C(i) + dt*incr;
    c(i+1) = C(i+1)/phi;
    
end

CDisc = C;
%% analytic solution
FNorm = F/hd;

%setup impuls-response function I
I = FNorm*exp(-FNorm/phi*timeline);

%setup convolution
C    = conv(I,aif)*dt;
CAna = C(1:k);



%% test deconvolution

%downsampling
step        = round(res/dt);
timelineMod = timeline(1:step:end);
CMod        = CDisc(1:step:end);
aifMod      = aif(1:step:end);
IMod        = I(1:step:end);
dtMod       = timelineMod(2)-timelineMod(1);

%prepare deconvolution
A       = perfusion1c.getLinearConvolutionMatrix(aifMod,dtMod);
[U,S,V] = svd(A);

%% do deconvolution
[FRec,IRec,CRec] = perfusion1c.linearDeconvolution(CMod,timelineMod,OI,U,S,V);

% show results

FTrue = FNorm;

if showDiscVsAna
    figure(1);clf;
    set(1,'name','Discrete vs. Analytic');
    plot(timeline,CAna,timeline,CDisc);
    legend('C_{Ana}','C_{Disc}');
end

if showDeconv
    figure(2);clf;
    set(2,'name','Deconvolution Results');
    
    subplot(1,2,1);
    plot(timelineMod,CMod,timelineMod,CRec);
    legend('CDisc','CRec');    
    title('C');
    
    subplot(1,2,2);
    plot(timelineMod,IMod,timelineMod,IRec);
    legend('I','IRec');
    ti = sprintf('FRec=%1.4f, FTrue=%1.4f',FRec,FTrue);
    title(ti);

end


