function [pmat,flowmat] = syntTPFALap(Kmat,Qmat,prm)
%function [pmat,flowmat] = syntTPFALap(Kmat,Qmat,prm)
% 
% Get the Pressure pmat by solving -K/mu * \Delta u = Q
%
% IT IS ASSUMED, THAT K IS CONSTANT IN OMEGA!!!
%
%
%                                          (c)Constantin Heck, 07-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
 
    m  = prm.dim(1:2);
    K  = max(Kmat(:));
    h  = prm.h(1:2);
    mu = prm.mu;


    %Dirichlet-BC on point (1,m(1)) (upper right corner to preserver symmetry)
    k  = m(1);


    %% setup laplacian matrix

    %function handle for identity
    id = @(i) speye(m(i));

    %setup finite difference matrix for laplacian in 1D
    Li = @(i) spdiags(ones(m(1)+1,1)*[1,-2,1],[-1,0,1],m(1),m(1)); 


    %get 1D-derivative-matrices with Neumann BC
    L1          = Li(1);
    L1(1,1)     = -1;
    L1(end,end) = -1;
    L1          = L1./h(1)^2;

    L2          = Li(2);
    L2(1,1)     = -1;
    L2(end,end) = -1;
    L2          = L2./h(2)^2;

    %kronecker to get the full Laplacian
    Lyy = kron(id(2),L1);
    Lxx = kron(L2,id(1));

    %Dirichlet boundary conditions on point (1,k)
    idxi = (k-1)*m(1)+1;
    idxj = (k-1)*m(1)+1;
    Lyy(idxi,idxj) = -2/h(1)^2;

    %setup Laplacian
    LAP = Lxx + Lyy;

    %% setup gradient (needed since q = \nabla p)

    %setup function handle for the gradient
    Di = @(i) spdiags(ones(m(1)+1,1)*[-1,1],[-1,0],m(1)+1,m(1)); 

    %get 1D-derivative matrices with Neumann BC
    D1          = Di(1);
    D1(1,1)     = 0;
    D1(end,end) = 0;
    D1          = D1/h(1);

    D2          = Di(2);
    D2(1,1)     = 0;
    D2(end,end) = 0;
    D2          = D2/h(2);

    %kronecker to get the full Laplacian
    GRADy = kron(id(2),D1);
    GRADx = kron(D2,id(1));


    %Dirichlet boundary conditions on point (1,k)
    idxi = (k-1)*m(1)+1;
    idxj = (k-1)*m(1)+1;
    GRADy(idxi,idxj) = 2/h(1);

    %% solve the system

    %setup operator and rhs
    A = -K/mu*LAP;
    b = Qmat;

    %solve the system
    p = A\b(:);

    %shift p, negative pressure is not physically plausible
    p = p-min(p(:));
    
    %reshape pmat
    pmat = reshape(p,m);

    %get the pressure as gradient
    flowmat{1} = reshape(GRADy*p(:),m(1)+1,m(2));
    flowmat{2} = reshape(GRADx*p(:),m(1),m(2)+1);
    flowmat{3} = zeros(m(1),m(2),1+1);

    
    
 
 
end 