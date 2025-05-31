function [ FRF, FBB ] = OMP( Fopt,N,M,AT,K)
    FRF_temp = complex(zeros(M*N,N));
    FBB_temp = complex(zeros(N,K));
    m = 1;
    Fres = Fopt;
    while m <= N(1) 
        Ns = m;
        Psi = AT'*Fres;
        [~,ind] = max(diag(Psi*Psi'));
        FRF_temp(:,m) = AT(:,ind);
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