function [F,Irec,Crec] = circularDeconvolution(C,AIF,timeline,OI,U,S,V) 
%function [F,Irec,Crec] = circularDeconvolution(C,AIF,timeline,OI) 
% 
% Circular Deconvolution according to [1].
% The convolution is modeled like in the function linearDeconvolution as:
%
%   C = M(AIF)*I
%
% The difference to linearDeconvolution is that 
% circularDeconvolution uses a circular convolution matrix M(AIF), i.e. 
% assumes periodic boundary conditions on AIF and on C. 
% The deconvolution is hence robust agains time-shifts of the AIF. In order
% to avoid aliasing-effects, AIF and C are zero-padded leading to a
% solution I twice the size of the original timeline.
% 
%
%
% INPUT:
%        C - Uptake Curve, size(C) = [k,1].
%      AIF - Arterial Input Function, size(AIF) = [k,1].
% timeline - timeline with measurments, size(timeline) = [k,1].
%
% VARARGIN:
%       OI - Threshold for Oszillation index. Default: OI=Inf.
%            See  also the documentation of linearDeconvolution.
%
% OUTPUT:
%        F - Flow, size(F)=[1,1];
%     Irec - F*R, size(Irec) = [2*k,1]
%     Cest - approximated uptake curve, size(Cest) = [2*k,1];
%     
%
% 
% REFERERENCES:
% [1] Wu, Ostergaard, Weisskoff et al.,
%     "Tracer arrival timing-insensitive technique for estimating flow in 
%      MR perfusion-weighted imaging using singular value decomposition 
%      with a block-circulant deconvolution matrix," 
%      MAGN RESON MED 50 (1),pp. 164-174, Jun. 2003.
%
%                                          (c)Constantin Heck, 24-Jan-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

%no input given? no problem!
if nargin == 0 
	clc
	help(mfilename) 
	runMinimalExample 
	return; 
end

persistent AIFsav;

%prepare variables
timeline = timeline(:);
AIF      = AIF(:);

%setup main variables
k        = numel(timeline);
CHat     = [C;zeros(k,1)];   %zero-padded C


AIFsav = AIF;


%prepare inverse of S
d       = diag(S);
dinv    = 1./d; 

%cut-off starts at first value where s~=0.
i = find(d<=eps,1);
if isempty(i), 
    i = 2*k+1; 
end;

%cut off until the current oscillation index OIc <= OI
OIc = Inf;
while OIc>OI
    
    %cut off small singular values
    dinv(i:end) = 0;
    
    %get the pseudoinverse of A
    Sinv  = diag(dinv,0);

    %recover Irec
    Irec  = V*(Sinv*(U'*CHat));
    
    %get oszillation index
    OIc = getOszillationIndex(Irec);
    i   = i-1;
    
end
    
    %get cutOff version of S
    cutOffd            = d;
    cutOffd((i+1):end) = 0; %i was already deminished in the loop
    cutOffS            = diag(cutOffd,0);
    
    %reconstruct
    Crec = U*(cutOffS*(V'*Irec));
    F    = max(Irec);
 

end 



function osc = getOszillationIndex(I)

    persistent H;
    n = numel(I);
    
    %matrix for oszillaton mearsurment
    if size(H,2)~=n
        H   = spdiags(ones(n,1)*[1,-2,1],[0,1,2],n-2,n); 
    end
    
    %get oszillation index
    osc = 1/(n*max(I))*sum(abs(H*I));

end


function runMinimalExample
    
    clear all;

    %setup the timeline and an AIF
    k        = 100;
    timeline = linspace(0,1,k)';

    %setup parameters for residue function = F*e^(-a*t)
    tau    = -.3; %parameter for shift
    MTT    = .4;
    F      = 20;
    deltaT = timeline(2)-timeline(1);


    %setup functions I and AIF
    AIF    = getGammaAIF(timeline,'t0',.4);
    AIFtau = getGammaAIF(timeline-tau,'t0',.4);
    I      = F*exp(-1/MTT*timeline);


    %forward model
    %setup matrix for forward model

    %get C
    C = conv(AIFtau,I);
    C = deltaT*C(1:k);

    
    %deconvolution
    OI = 3e-2; %determined experimentally
    [Frec,Irec,Crec] = circularDeconvolution(C,AIF,timeline,OI);

    
    %plot results
    %setup extended timeline timeline2
    timeline2 = linspace(0,2,2*k);


    %plot results
    figure(1);clf;
    subplot(1,3,1);
    plot(timeline,AIF,timeline,AIFtau)
    legend('AIF','delayed AIF');

    subplot(1,3,2);
    plot(timeline,AIFtau,timeline,C,timeline,I,timeline2,Crec);
    legend('delayed AIF','C','I','Crec')
    title(sprintf('Error in Flow: |F-Frec|/F=%1.2e',abs(F-Frec)/F))
    

    subplot(1,3,3);
    plot(timeline,I,timeline2,Irec);
    legend('I','Irec')
    
    F
    Frec


end