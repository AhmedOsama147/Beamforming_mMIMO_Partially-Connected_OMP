function [FRF,FBB] = SIC(K,N,M,NR,H,SNR)
        T0 = eye(NR);
        R  = [eye(M) zeros(M,M*(N-1))];
        G0 = R*H'*inv(T0)*H*R';
        G = G0; 
        A = zeros(M*N,N); D = zeros(N,N);
        EN = 10.^(SNR/10);
        for m = 1:N
            [~,S,V] = svd(G);
            v1 = V(:,1); s = S(1,1);
            ag = (1/sqrt(M))*norm(v1,1)*exp(1i*angle(v1));
            A((m-1)*M+1:m*M,m) = ag;
            dg = (1/sqrt(M))*norm(v1,1);
            D(m,m) = dg;
            X = ((EN/K)*(s^2)*v1*v1')/(1+(EN/K)*s);
            G = G - X;            
        end
        FRF = A ; FBB = D(:,1:K);
        Fnf = sqrt(trace((FRF*FBB)*(FRF*FBB)'));
        FBB = sqrt(K)*FBB/Fnf;
end