clc;clear; close all

Nray = 10;
Ncl = 8;
maxCap = 22;
CapStep = 2;
minCap = 6;
minNR = 8;
maxNR = 32;
steps = 2;
NR = minNR:steps:maxNR;
SNR = 0;
EN = 10.^(SNR/10);
lenNR = length(NR);
  
NT_RF  = 8; %RF chain at TX
NT = 64; %Antenna elements at TX
NS = 3; %Data Streams

num = 100;

alpha = sqrt(1/(2)) *  (randn(Ncl,Nray,num) + 1i * randn(Ncl,Nray,num)) ;
phiT = (-pi) + ((2*pi)*rand(Ncl,num));
phiR = (-pi) + ((2*pi)*rand(Ncl,num));
sig = deg2rad(7.5); % this is the s.d. of the Laplace distribution for the scattering clusters

Capacity_fully_digital = zeros(num,1);
cap_full_digital = zeros(lenNR,1);
Capacity_OMP = zeros(num,1);
cap_OMP = zeros(lenNR,1);
Capacity_NKP = zeros(num,1);
cap_NKP = zeros(lenNR,1);
Capacity_SIC = zeros(num,1);
cap_SIC = zeros(lenNR,1);
Capacity_mod = zeros(num,1);
cap_mod = zeros(lenNR,1);

for k = 1:lenNR
    ff = zeros(num,1);
    for i = 1:num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        At = [];
        Phi_AOD = zeros(Ncl,Nray);
        Phi_AOA = zeros(Ncl,Nray);
        H = zeros(NR(k),NT);
        for j=1:Ncl
             % The final azimuth angles from each cluster i    
             phi_AOD = randLaplacian(Nray, 1, phiT(j,i), sig^2);
             phi_AOA = randLaplacian(Nray, 1, phiR(j,i), sig^2);    
             for l=1:Nray
                 at = exp(-1j*pi*sin(phi_AOD(l))*(0:1:NT-1)).'/sqrt(NT);
                 ar = exp(-1j*pi*sin(phi_AOA(l))*(0:1:NR(k)-1)).'/sqrt(NR(k));
                 Alpha = alpha(j,l,i);
                 H = H + Alpha*ar*at'; % channel
                 At = [At at];
             end
         end
         H_n = sqrt(NT*NR(k)/(Ncl*Nray))*H;
         ff(i) = norm(H_n,'fro')^2;
         [U,~,V] = svd(H_n);
         Fopt = V(:,1:NS);
         Capacity_fully_digital(i) = log2(det(eye(NR(k))+(EN/NS)*H_n*Fopt*Fopt'*H_n'));
         
         [ FRF, FBB ] = OMP(Fopt,NT_RF,NT/NT_RF,At,NS);
         Capacity_OMP(i) = log2(det(eye(NR(k))+(EN/NS)*H_n*FRF*FBB*FBB'*FRF'*H_n'));
        
         [FRF_NKP,FBB_NKP] = NKP(NS,NT_RF,NT/NT_RF,Nray,Fopt);
         Capacity_NKP(i) = log2(det(eye(NR(k))+(EN/NS)*H_n*FRF_NKP*FBB_NKP*FBB_NKP'*FRF_NKP'*H_n'));
         
         [FRF_SIC,FBB_SIC] = SIC(NS,NT_RF,NT/NT_RF,NR(k),H_n,SNR);
         Capacity_SIC(i) = log2(det(eye(NR(k))+(EN/NS)*H_n*FRF_SIC*FBB_SIC*FBB_SIC'*FRF_SIC'*H_n'));
         
         [FRF_mod,FBB_mod] = Mod_OMP(NS,NT_RF,NT/NT_RF,Fopt,At);
         Capacity_mod(i) = log2(det(eye(NR(k))+(EN/NS)*H_n*FRF_mod*FBB_mod*FBB_mod'*FRF_mod'*H_n'));
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cap_full_digital(k) = mean(abs(Capacity_fully_digital));
cap_OMP(k) = mean(abs(Capacity_OMP));
cap_NKP(k) = mean(abs(Capacity_NKP));
cap_SIC(k) = mean(abs(Capacity_SIC));
cap_mod(k) = mean(abs(Capacity_mod));


end

plot(NR,smooth(cap_full_digital),'-ks','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','black','MarkerFaceColor',[.0 .0 .0]);
hold on;
plot(NR,smooth(cap_OMP),'--ms','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','m','MarkerFaceColor','m');
hold on;
plot(NR,smooth(cap_NKP),'-rs','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
hold on;
plot(NR,smooth(cap_SIC),'-bs','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','blue','MarkerFaceColor',[.3 .5 .8])
hold on;
plot(NR,smooth(cap_mod),'-gs','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','green','MarkerFaceColor',[0 .5 .5]);
hold on;
xlim([minNR maxNR]);
ylim([minCap maxCap]);
grid on;
% title("Spectral efficiencies achieved by the proposed hybrid analog-digital precoding algorithm in a 64×16 mmWave system");
xticks(minNR:steps:maxNR);
yticks(minCap:CapStep:maxCap);
set(get(gca,'XLabel'),'String','R_{MS}', 'fontname', 'Times');
set(get(gca,'YLabel'),'String','Capacity (bps/Hz)','Interpreter','latex');
grid on;

h1 = legend('Fully digital Precoding',...
            'OMP-Fully connected','NKP-Partially connected','SIC-Partially connected',...
            'Proposed Algorithm',...
           'Location','Northwest');
set(h1, 'Fontsize', 8,'Interpreter','latex');