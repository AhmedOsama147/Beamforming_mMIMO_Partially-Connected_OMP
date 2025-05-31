function [FRF,FBB] = NKP(K,N,M,L,Fopt)
    F = [];
    for j = 1:K
        Q = Fopt(:,j);
        Q1 = reshape(Q,[M N]);
        Q1 = transpose(Q1);
        F = [F;Q1];
    end
    [~,~,V] = svd(F);
    v1 = V(:,1);
    angV = angle(v1);
    f_RF = (1/sqrt(M))*exp(-1i*angV);
    FRF = kron(eye(N),f_RF);
    FBB = (FRF'*Fopt);
    Fnf = sqrt(trace((FBB)*(FBB)'));
    FBB = sqrt(K)*FBB/Fnf;    
end