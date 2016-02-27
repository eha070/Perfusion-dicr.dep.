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
% Then recovers F by deconvolving F.
% 
%                                      (c)Constantin Sandmann, 27-Feb-2016 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


showDiscVsAna = 1;
showDeconv    = 1;

%% setup main variables

%setup main constants
F   = 50/60; %flow in ml^3/s/ml
hd  = 100;   %voxel volume in ml
phi = 0.05;  %porosity in ml/ml


%main input parameters
timeline = linspace(0,90,1001)';
aif      = perfusion1c.getGammaAIF(timeline/60);
dt       = timeline(2)-timeline(1);
k        = numel(timeline);


%oscillation index
OI = .003;

%% solution by discretization

c = zeros(k,1);
C = zeros(k,1);
for i = 1:k-1
    
    incr = F*(aif(i) - c(i));
    
    C(i+1) = C(i) + dt*(incr/hd);
    c(i+1) = C(i+1)/phi;
    
end

CDisc = C;
%% analytic solution
F_Norm = F/hd;

%setup impuls-response function I
I = F_Norm*exp(-F_Norm/phi*timeline);

%setup convolution
C     = conv(I,aif)*dt;
CAna = C(1:k);



%% test deconvolution

%downsampling
step        = round(1/dt);
timelineMod = timeline(1:step:end);
CMod        = CDisc(1:step:end);
aifMod      = aif(1:step:end);
IMod        = I(1:step:end);
dtMod       = timelineMod(2)-timelineMod(1);

%prepare deconvolution
A       = perfusion1c.getLinearConvolutionMatrix(aifMod,dtMod);
[U,S,V] = svd(A);

%do deconvolution
[FRec,IRec,CRec] = perfusion1c.linearDeconvolution(CMod,timelineMod,OI,U,S,V);

%setup true F
FTrue = F;


%% show results

if showDiscVsAna
    figure(1);clf;
    set(1,'name','Discrete vs. Analytic');
    plot(timeline,CAna,timeline,CDisc);
    legend('C_{Ana}','C_{Disc}');
end

if showDeconv
    figure(2);clf;
    set(1,'name','Deconvolution Results');
    
    subplot(1,2,1);
    plot(timelineMod,CMod,timelineMod,CRec);
    legend('C_Disc','CRec');    
    title('C');
    
    plot(timelineMod,IMod,timelineMod,IRec);
    legend('I','IRec');
    ti = sprintf('F=%1.4f, FTr=%1.4f',F,FTrue);
    title(ti);

end


