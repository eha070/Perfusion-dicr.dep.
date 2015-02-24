

function A = getLinearConvolutionMatrix(AIF,deltaT)
    
    %setup main variables
    AIF = AIF(:);
    n = numel(AIF);

    %standard blockcirculant matrix
    A  = AIF;
    for i = 1:(n-1)
        A      = [A,circshift(AIF,i,1)];
    end
    
    %make it a lower-triangular matrix
    A = tril(A,0);
    
    %include integration
    A = deltaT*A;



end
