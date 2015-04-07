function F = maximumSlope(Cmat,timeline,aif) 
%function F = maximumSlope(Cmat,timeline,aif)  
% Uses the maximum slope technique to estimate CBF as implemented for 
% SIEMENS scanners. Fick's law is assumed and CBF is estimated as
%
%   CBF = max(d/dt C)/max(aif).
% 
% INPUT:
%      Cmat - the tissue curves.
% timeline  - timeline for measurments
%      aif  - the Arterial Input Function
%
% OUTPUT:
%      F  - The map for CBF
% 
% REFERENCES:
% [1] Klotz, Koenig, "Perfusion measurements of the brain: using dynamic CT 
%     for the quantitative assessment of cerebral ischemia in acute stroke"
%     1999, European Journal of Radiology(30),170-?184
%
%                                          (c)Constantin Heck, 24-Mar-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
 
    if nargin == 0 
        clc
        help(mfilename) 
        return; 
    end

    %initialize variables
    k        = numel(aif);
    aif      = aif(:);
    Cmat     = reshape(Cmat,[],k);
    Cmat     = Cmat';
    timeline = timeline(:);

    %setup matrix for derivatives
    tau = 1./(timeline(2:end)-timeline(1:end-1));
    D = spdiags(tau*[-1,1],[0,1],k-1,k);
      
    %get maximum slope
    MS   = max(D*Cmat,[],1);
    maif = max(aif);
    
    %get flow map
    F = MS./maif;
    

 
end 