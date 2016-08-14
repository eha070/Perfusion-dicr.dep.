function varargout = getParkerAIF(timeline,varargin)
%function varargout = ParkerAIF(timeline,varargin)
%   INPUT:
%   timeline - timepoints to evaluate the AIF (min).
%
%   VARARGIN:
%   p           - parametervector, p \in \R^{10}. Default: p from Parker's paper.
%   dAIFdt   - 0 or 1. Default: 0. Returns derivative w.r.t. t, i.e.
%               - varargout = {AIF,dAIF,dAIFdt}
%
%   OUTPUT:
%   AIF - The AIF in [mmol/l]. Concentration in the blood.
%   dAIF - the derivative of the AIF w.r.t. p. dAIF has dimension length(t) \times 10
%   dAIFdt - the derivative of the AIF w.r.t. t. dAIF has dimensions
%            length(t) \times length(t)


%display help and run minimal example
if nargin == 0 
  help(mfilename);
  disp('-------------------')
  disp('Default parameters by Parker are used')
  p = [ 0.809, 0.330, 0.17046, 0.365, 0.0563, 0.132, 1.050, 0.1685, 38.078, 0.483]';
  timeline = (0:.01:6)';
  AIF = getParkerAIF(timeline,'p',p);
  plot(timeline,AIF);
  return;
end

numMeas = length(timeline);
timeline = timeline(:);
p = [ 0.809, 0.330, 0.17046, 0.365, 0.0563, 0.132, 1.050, 0.1685, 38.078, 0.483]';
dAIFdt = 0;

% overwrite defaults 
for k=1:2:length(varargin)
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

p = reshape(p,[],1);


A1 = p(1);
A2 = p(2);
T1 = p(3);
T2 = p(4);
sigma1 = p(5);
sigma2 = p(6);
alpha = p(7);
beta = p(8);
s = p(9);
tau = p(10);






C = @(t) A1/(sigma1*sqrt(2*pi))*exp( -(t-T1).^2/(2*sigma1^2) ) ...
       + A2/(sigma2*sqrt(2*pi))*exp( -(t-T2).^2/(2*sigma2^2) ) ...
       + alpha*exp(-beta*t)./( 1 + exp( -s*(t - tau) ) );

dCdA1       = @(t) 1/(sigma1*sqrt(2*pi))*exp( -(t-T1).^2/(2*sigma1^2) );
dCdA2       = @(t) 1/(sigma2*sqrt(2*pi))*exp( -(t-T2).^2/(2*sigma2^2) );
dCdT1       = @(t) A1/(sigma1*sqrt(2*pi))*2/(2*sigma1^2)*(t-T1).*exp( -(t-T1).^2/(2*sigma1^2) );
dCdT2       = @(t) A2/(sigma2*sqrt(2*pi))*2/(2*sigma2^2)*(t-T2).*exp( -(t-T2).^2/(2*sigma2^2) );
dCdSigma1   = @(t) -A1/(sigma1^2*sqrt(2*pi))*exp( -(t-T1).^2/(2*sigma1^2) ) + A1/(sigma1^4*sqrt(2*pi))*(t-T1).^2.*exp( -(t-T1).^2/(2*sigma1^2) );
dCdSigma2   = @(t) -A2/(sigma2^2*sqrt(2*pi))*exp( -(t-T2).^2/(2*sigma2^2) ) + A2/(sigma2^4*sqrt(2*pi))*(t-T2).^2.*exp( -(t-T2).^2/(2*sigma2^2) );
dCdAlpha    = @(t) exp(-beta*t)./( 1 + exp( -s*(t - tau) ) );
dCdBeta     = @(t) alpha*(-t).*exp(-beta*t)./( 1 + exp( -s*(t - tau) ) );
dCdS        = @(t) -alpha*exp(-beta*t).*( 1 + exp( -s*(t - tau) ) ).^(-2).*( -(t-tau) ).*exp( -s*(t - tau) );
dCdTau      = @(t) -alpha*exp(-beta*t).*( 1 + exp( -s*(t - tau) ) ).^(-2).*exp( -s*(t - tau) )*s;



AIF =   C(timeline);
dAIF =  [ dCdA1(timeline), dCdA2(timeline), dCdT1(timeline), dCdT2(timeline), dCdSigma1(timeline), dCdSigma2(timeline), dCdAlpha(timeline), dCdBeta(timeline), dCdS(timeline), dCdTau(timeline) ];


if dAIFdt == 0
    varargout = {AIF,dAIF};
    
elseif dAIFdt== 1
    dCdt = @(t) A1/(sigma1*sqrt(2*pi))*exp( -(t-T1).^2/(2*sigma1^2) ).*(-(t-T1)/sigma1^2) ...
       + A2/(sigma2*sqrt(2*pi))*exp( -(t-T2).^2/(2*sigma2^2) ).*(-(t-T2)/sigma2^2) ...
       - alpha*exp(-beta*t)./(( 1 + exp( -s*(t - tau) ) ).^2).*( -s*exp( -s*(t-tau) ) ) ...
       + alpha*(-beta)*exp(-beta*t)./( 1 + exp( -s*(t - tau) ) );
    varargout = {AIF,dAIF,diag(dCdt(timeline))};
    
end


end