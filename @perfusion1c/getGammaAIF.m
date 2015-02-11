function AIF = getGammaAIF(timeline,varargin) 
% function AIF = getGammaAIF(timeline,varargin) 
% Sets up an AIF according to the model
% 
%     AIF(t) = C0*(t-t0)^a*exp(-(t-t0)/b) for t>=t0
%     AIF(t) = 0                          for t<t0
% with a=3.0 and b=1.5s. The timeline is supposed to be given in minutes.
%
% INPUT:
%   timeline - Given in minutes
%
% VARARGIN:
%         t0 -   Delay (in min): default 0.0
%         C0 - scaling         : default 1.0
%
% OUTPUT:
%        AIF - The AIF.
% 
%   
%
% REFERENCES:
% [1] Wu, Ostergaard, Weisskoff et al.,
%     "Tracer Arrival Timing-Insensitive Technique for Estimating Flow in 
%     MR Perfusion-Weighted Imaging Using Singular Value Decomposition With 
%     a Block-Circulant Deconvolution Matrix",
%     MAGN RESON MED 50, (p.164--174), 2003 and references within.
%
%                                          (c)Constantin Heck, 05-Dec-2014 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
%if no input is given: Show the AIF.
if nargin == 0 
	clc
	help(mfilename) 
	runMinimalExample 
	return; 
end

%default parameters
C0 = 1;
t0 = 0;
 
%overwrites default parameter 
for j=1:2:length(varargin),  
	eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
end 

%fixed parameters
a  = 3.0;
b  = 1.5;


%change from minutes to seconds
t  = timeline*60;
t0 = t0*60;


%get the AIF
AIF = C0*(t-t0).^a.*exp(-(t-t0)./b);
AIF = AIF.*(t>=t0);
 
 
end 
 
 
function runMinimalExample 

    %get timeline in minutes
    timeline = (0:.001:1);
    t0       = 0.0; %delay in minutes
 
    %get the AIF
	gammaAIF = feval(mfilename,timeline,'t0',t0);
    
    %plot result
    figure(1);clf;
    plot(timeline,gammaAIF);
    legend('gammaAIF');
    xlabel('time (min)');
    ylabel('concentration');
    title(sprintf('Gamma-variate AIF from Wu et al. (2003),\n delay=%1.2fmin',t0));
end 