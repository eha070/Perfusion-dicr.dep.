function [F,Irec,Crec] = circularDeconvolution(C,timeline,OI,U,S,V) 
%function [F,Irec,Crec] = circularDeconvolution(C,timeline,OI,U,S,V) 
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
        return; 
    end

    %prepare variables
    timeline = timeline(:);

    %setup main variables
    k        = numel(timeline);
    CHat     = [C;zeros(k,1)];   %zero-padded C


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
