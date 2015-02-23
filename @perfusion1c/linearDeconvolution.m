function [F,Irec,Crec] = linearDeconvolution(C,timeline,OI,U,S,V) 
%function [F,Irec,Crec] = linearDeconvolution(C,timeline,OI,U,S,V) 
% NEW VERSION
% 
% Linear Deconvolution according to the model:
%
%   C = M(AIF)*I
%
% where M(AIF) is a lower-triangular matrix. I is then recovered by
%   
%   I = M^\dagger * C
% 
% where M^\dagger is the pseudo-inverse of M, obtained by a thresholded
% SVD. The general technique is described in [1].
% The threshold dependends on the value of OI, the oscillation index. 
%
%   O(I) = scaled-version of the second derivatives of I
%
% the threshold for the SVD is increased until O(I)<OI.
%
%
%
% INPUT:
%        C - Uptake Curve
%      AIF - Arterial Input Function
% timeline - timeline with measurments
%
% VARARGIN:
%       OI - Threshold for Oszillation index. Default: OI=Inf;
%
% OUTPUT:
%        F - Flow
%        R - Residue function
%     Irec - F*R
%     Cest - approximated uptake curve
%     
%
% REFERENCES:
% [1] Ostergaard, Weiskoff, Chesler, Gyldensted, Rose,
%     "High resolution measurement of cerebral blood flow using 
%      intravascular tracer bolus passages. Part I: Mathematical approach 
%      and statistical analysis", 
%      Magn Reson Med 36(5), p.715-725, 1996
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



    %prepare inverse of S
    d       = diag(S);
    dinv    = 1./d; 

    %cut-off start from value i
    i = find(d<=eps,1);
    if isempty(i), 
        i = 2*k+1; 
    end;

    %cut off until osc < OI
    osc = Inf;
    while osc>OI

        %cut off small singular values
        dinv(i:k) = 0;

        %get the pseudoinverse of A
        Sinv  = diag(dinv,0);

        %recover Irec
        Irec  = V*(Sinv*(U'*C));

        %get oszillation index
        osc = getOszillationIndex(Irec);
        i   = i-1;

    end

    %get cutOff version of S
    cutOffd            = d;
    cutOffd((i+1):end) = 0; %i was already deminished in the loop
    cutOffS            = diag(cutOffd,0);

    
    %get variables
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
    osc = 1/(2*n*max(I))*sum(abs(H*I));

end