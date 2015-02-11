function [name] = providenameindicator(phiopt,Kopt,dim,aiftype,T)

str = perfusion1c.providenameflow(phiopt,Kopt,dim);
name = [str '-aif-' aiftype '-T-' num2str(T)];

