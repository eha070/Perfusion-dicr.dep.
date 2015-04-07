classdef perfusion1c
       
   methods (Static)
      % main functions
      createflowTPFA;
      [prm,Qmat] = settings;
      [pmat,qmat] = syntTPFA(Kmat,Qmat,prm);
      [perfmat,perfmatn] = flux2perf(qmat,lenim,Qmat,h);
      basenameflow = providenameflow(phiopt,Kopt,dim);
      basenameindicator = providenameindicator(phiopt,Kopt,dim,aiftype,T);
      basenamerecon = providenamerecon(phiopt,Kopt,dim,aiftype,T,deconvmethod,red);
      createindicatorpde(prm);
      Cmat = syntforwprobpde(phimat,qmat,Qmat,aifval,prm);
      aifval = createaif(aiftype,timeline);
      reducetimesteps(pathload,stepred);
      createindicatorconv(prmin);
      createindicatorconvav(prmin);
      [qmat,perfmat,perfmatn,phimat] = reconflow(prmin);
      [phimat] = porosity(Cmat,aifval,timeline);
      [qmat,perfmat,perfmatn] = reconflowclassic(Cmat,timeline,roi,aifval,prm);
      aifval = getGammaAIF(timeline);
      aifval = getParkerAIF(timeline);
      A = getCircularConvolutionMatrix(AIFHat,deltaT);
      A = getLinearConvolutionMatrix(AIF,deltaT);
      [F,Irec,Crec] = circularDeconvolution(C,timeline,OI,U,S,V); 
      [F,Irec,Crec] = linearDeconvolution(C,timeline,OI,U,S,V);
      [lenim] = arclength(q,Q,h);
      [a] = distimage(minima,conn);
      [im] = constborder(im,p,v);
      [v] = interp3c(im,c);
      [prmdefault] = mergestruct(prmdefault,prmin);
      [h] = panelstruct(varargin)
      printstructscreen(var);
      [data] = im2vec4D(varargin);
      [A] = mat2celldirect(a);
      cell2tex(A,filename,n);
      [u] = transim(varargin);
      basename = struct2nameIndicator(prm,varargin);
      qcc = convertFlowStagToCC(qmat);
      F = maximumSlope(C,timeline,aif);
      CBV = cbvEstimation(C,timeline,aif);

   end 
      
end
