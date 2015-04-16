function pOpt = FitGaussian(timeline,C,p0,varargin) 
%function pOpt = FitGaussian(timeline,C,p0,varargin)  
% 
% Fits a 1D-Gaussian function to timeline,C with a starting guess p0
%
%   M = a*exp(-((x - b)/c)^2)
%
%
% 
% INPUT:
% timeline - the timlein
%        C - the concentration curve
%       p0 - starting guess
%
% VARARGIN:
%   maxIter,tolG,tolY,tolJ - Optimization paramters
%                 showIni  - only show the initial fit
%
% OUTPUT:
%     pOpt = Optimized parameters
%      
% 
%                                          (c)Constantin Heck, 12-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
 
%what to do if no input is given
if nargin == 0 
	clc
	help(mfilename) 
	runMinimalExample 
	return; 
end

%default parameters
maxIter = 30;
tolG    = 1e-3;
tolY    = 1e-5;
tolJ    = 1e-3;
showIni = false;

 

%overwrites default parameter 
for j=1:2:length(varargin),  
	eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
end 


%setup objective function and get initial model-fit
fun  = @(p) obFun(p,timeline,C);
M0   = getModel(p0,timeline);

%if showIni: Show initial guess,quit. Else.... well?
if showIni
    figure(1);clf;
    plot(timeline,C,timeline,M0);
    legend('C','M0');
    
    pOpt = [];
    return;
    
else
    %optimization
    pOpt = GaussNewtonCheck(fun,p0,'maxIter',maxIter,'tolG',tolG,'tolY',tolY,'tolJ',tolJ);
    
    %get results
    [~,para] = fun(pOpt);
    M        = para.M;

    %show results
    figure(1);clf;
    plot(timeline,C,timeline,M);
    legend('C','M');    
end









 
 
 
 
 
 
end 



function [f,para,df,d2f] = obFun(p,timeline,C)

    [M,dM] = getModel(p,timeline);
    

    res = M - C;    
    f   = 1/2*(res'*res);
    df  = res'*dM;
    d2f = dM'*dM;
    
    para = struct('M',M);
    


end


 
function [M,dM] = getModel(p,timeline)
%function [M,dM] = getModel(p,timeline)
%
%Simple 1D-Gauss-Distr:
% 
% M = ae^(- ((x - b)/c)^2)
%


    a = p(1);
    b = p(2);
    c = p(3);
    
    timeline = timeline(:);
    
    
    

    E = exp(-((timeline - b)./c).^2);
    e = -2*(timeline - b)./c;
    M = a*E;    
    
    
    dMda = E;
    dMdb = M.*e.*(-1/c);
    dMdc = M.*e.*(-(timeline-b)./c^2);
    dM   = [dMda,dMdb,dMdc];


end

 
function runMinimalExample 

    timeline = linspace(-5,5,100);
    pIni     = rand(3,1);
    p0       = rand(3,1);
    C        = getModel(pIni,timeline);

%     %check model
%     g  = @(p) getModel(p,timeline);
%     checkDerivative(g,p0);

%     %check L2-distance
%     g  = @(p) obFun(p,timeline,C);
%     checkDerivative(g,p0);
    
    
    %check the whole function
    p0   = pIni + .05*randn(size(pIni));
    pOpt = feval(mfilename,timeline,C,p0);
    
    %does it work?
    pIni - pOpt
    

end 