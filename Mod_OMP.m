function [FRF,FBB] = Mod_OMP(K,N,M,Fopt,AT)
    FRF_temp = complex(zeros(M*N,N));
    FBB_temp = complex(zeros(N,K));
    m = 1;
    Fres = Fopt;
    while m <= N(1) 
        Ns = m;
        A = AT((m-1)*M+1:m*M,:);
        Psi = A'*Fres((m-1)*M+1:m*M,:);
        [~,ind] = max(diag(Psi*Psi'));
        FRF_temp((m-1)*M+1:m*M,m) = (1/sqrt(M))*exp(1i*angle(A(:,ind)));
        FBB_temp(1:m,:) = pinv(FRF_temp(:,1:m)'*FRF_temp(:,1:m))*(FRF_temp(:,1:m)'*Fopt);
        temp = Fopt-FRF_temp(:,1:m)*FBB_temp(1:m,:);
        Fn = sqrt(trace((Fres - FRF_temp*FBB_temp)*(Fres - FRF_temp*FBB_temp)'));
        Fres = temp/Fn;
        m = m+1;
    end
    FRF = FRF_temp(:,1:Ns);
    FBB = FBB_temp(1:Ns,:);
    Fnf = sqrt(trace((FRF*FBB)*(FRF*FBB)'));
    FBB = sqrt(K)*FBB/Fnf;
end