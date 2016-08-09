function T = getLinearInterMatrix1D(x,xq,varargin) 
%function A = getLinearInterMatrix1D(x,xq,varargin)  
% Sets up a matrix T such that T*v=interp1(x,v,xq)
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 13-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
 
if nargin == 0 
	clc
	help(mfilename) 
	runMinimalExample 
	return; 
end
 
 
%overwrites default parameter 
for j=1:2:length(varargin),  
	eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
end 


 
 

%prepare variables
x  = x(:);
xq = xq(:);
nx = numel(x);
nq = numel(xq);






%{
%setup function
e =@(i) double((1:nx) == i)';

%setup matrix
T = zeros(nq,nx);
for i = 1:nx
    
    ei = e(i);
    f  = interp1(x,ei,xq,'linear');
    T(:,i) = f(:);
    
end

T(isnan(T)) = 0;
%}




%find vector I = (j_1,...,j_nq) such that x(j_i) <= xq(i)
I = zeros(nq,1);
for i = 1:nq
    xqi = xq(i);
    
    I(i) = nx - sum(xqi<x);

end


%find indicies which are not valid
val  = xq<max(x) & xq>=min(x);
nval = nnz(val);
Ival = I(val);

%this index is valid but I(val)+1 is out of range
valm = (xq==max(x));



%setup weights
z = (x(Ival+1) - xq(val))./(x(Ival+1)-x(Ival));




%setup indices
Tz = zeros(nval,nx);
idxz0 = sub2ind(size(Tz),(1:nval)',Ival);
Tz(idxz0) = z;
idxz1 = sub2ind(size(Tz),(1:nval)',Ival+1);
Tz(idxz1) = (1-z);

%setup matrix
T        = zeros(nq,nx);
T(val,:) = Tz;

%deal with out of range index
T(valm,end) = 1;
 
 
end 
 
 
function runMinimalExample 
	

    n1 = 5;
    n2 = 100;
    
    x  = linspace(0,1,n1);
    xq = linspace(0,2,n2);
    
    
    v  = rand(n1,1);
    T  = feval(mfilename,x,xq);
    
    Tv = T*v;
    
    
    figure(1);clf;
    plot(x,v,'-b',xq,Tv,'-r','LineWidth',2,'MarkerSize',5);
    legend('Given Data','Interpolated Data');
    
    inter = interp1(x,v,xq)';
    inter - Tv
 
end 