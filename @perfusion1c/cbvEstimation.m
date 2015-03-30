function CBV = cbvEstimation(C,timeline,aif) 


%setup main parameters
k        = numel(timeline);
timeline = timeline(:);
C        = C(:);
aif      = aif(:);

%get vector for integration
e  = timeline(2:end)-timeline(1:end-1);
Av = spdiags(1/2*ones(k,2),[0,1],k-1,k);
e  = e'*Av; e=e(:);


CBV = (e'*C)/(e'*aif(:));










end
