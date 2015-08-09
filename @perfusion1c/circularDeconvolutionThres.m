function [F,Irec,Crec] = circularDeconvolutionThres(C,timeline,thres,U,S,V) 
%function [F,Irec,Crec] = circularDeconvolutionThres(C,timeline,thres,U,S,V) 
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
% DIFFERENCE TO circularDeconvolution:
% Here the deconvolution is done with a fixed threshold for ALL curves.              
% 
%
%
% INPUT:
%        C - Uptake Curve, size(C) = [k,1].
%      AIF - Arterial Input Function, size(AIF) = [k,1].
% timeline - timeline with measurments, size(timeline) = [k,1].
%
% VARARGIN:
%    thres - Relative threshold for the eigenvalues w.r.t. the larges
%            eigenvalue.
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
    CHat     = [C(:);zeros(k,1)];   %zero-padded C


    %prepare inverse of S
    d       = diag(S);
    dinv    = 1./d; 
    
    %cut-off all eigenvalues which are too large
    thresEig  = thres*max(d);
    idx       = (d<=thresEig);
    
    %get Sinv
    dinv(idx) = 0;
    
    %get cutOffS
    d(idx)    = 0;


    %recover Irec
    Irec  = V*(dinv.*(U'*CHat));

    
    %reconstruct
    Crec = U*(d.*(V'*Irec));
    F    = max(Irec);

end 