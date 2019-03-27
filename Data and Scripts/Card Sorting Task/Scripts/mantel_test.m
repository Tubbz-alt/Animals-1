
function [zM,zN,rho,p] = mantel_test(A,B,nPerms)
% modified from http://richiardi.net/code/jMantelTest.m
% written by Judy Sein Kim 

% zM: Mantel test stat (cross-product of individual matrix entries aka sum of the products of corresponding elements) 
% zN: normalized Mantel test statistic (same as Pearson's r). 
% rho: Spearman's rho (we use this instead of zn because sorting data are very non-normal (large clusters around 0, problem with converting count data to %). 
% p: null hypothesis is no association between matrices ("distance between pairs indepedent of spatial distance between pairs")
% Notes: don't confuse zM with standard z-score (it is not standardized) 

    [I] = itril(size(A,1),-1); 
    Av = A(I); 
    Bv = B(I); 
    
    zM = Av'*Bv; % Reporting this for comparison, but not used 
    zN = corr(Av,Bv,'type','Pearson'); % Reporting this for comparison, but not used 
    rho = corr(Av,Bv,'type','Spearman');
    rhoPerm = zeros(nPerms,1); 
    
    for ind = 1:nPerms 
        % Permute rows and columns of B 
        randInd = randperm(size(B,1)); 
        BPerm = B(randInd,randInd); 
        BvPerm = BPerm(I);
        rhoPerm(ind) = corr(Av,BvPerm,'type','Spearman'); 
    end
    p = (1+sum(abs(rhoPerm)>=abs(rho)))/(nPerms+1); 
    figure
    hist(rhoPerm) 
    hold on; 
    line([rho,rho],ylim,'LineWidth',2,'Color','r'); 
end 
