function ret = Kapur(N,D,xR,Prob)

ret = zeros(1,N);

for j = 1: N
    PI0 = Prob(1:xR(j,1)); 
    ind = PI0 == 0;
    ind = ind .* eps;
    PI0 = PI0 + ind;
    
    w0 =  sum(PI0); 
    H0 = -sum((PI0/w0).*(log2(PI0/w0)));
    ret(j) = ret(j) + H0;
    
    for jl = 2: D
        PI = Prob(xR(j,jl-1)+1:xR(j,jl));
        ind = PI == 0;
        ind = ind .* eps;
        PI = PI + ind;
        
        w =  sum(PI); 
        H = -sum((PI/w).*(log2(PI/w)));
        ret(j) = ret(j) + H;
    end  
end

ret = ret.';
