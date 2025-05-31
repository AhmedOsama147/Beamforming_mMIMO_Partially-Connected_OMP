clc;clear; close all

Nray = 10;
Ncl = 8;
maxCap = 40;
CapStep = 5;
minCap = 0;
minSNR = -30;
maxSNR = 20;
steps = 5;
SNR = minSNR:steps:maxSNR;
EN = 10.^(SNR/10);
lenSNR = length(SNR);
  
NT_RF  = 2.^(0:6); %RF chain at TX
lenN = length(NT_RF);
NT = 64; %Antenna elements at TX
NR = 16; %Antenna elements at RX
NS = 3; %Data Streams

num = 100;
gamma = sqrt((NT*NR)/(Ncl*Nray));
alpha = sqrt(1/(2)) *  (randn(Ncl,Nray,num) + 1i * randn(Ncl,Nray,num)) ;
phiT = (-pi) + ((2*pi)*rand(Ncl,num));
phiR = (-pi) + ((2*pi)*rand(Ncl,num));
sig = deg2rad(7.5); % this is the s.d. of the Laplace distribution for the scattering clusters


Capacity_mod = zeros(num,lenSNR,lenN);
cap_mod = zeros(lenSNR,lenN);
Capacity_fully_digital = zeros(num,1);
cap_full_digital = zeros(lenSNR,1);

for k = 1:lenSNR
    for m = 1:lenN
    ff = zeros(num,1);
    for i = 1:num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        At = [];
        Phi_AOD = zeros(Ncl,Nray);
        Phi_AOA = zeros(Ncl,Nray);
        H = zeros(NR,NT);
        for j=1:Ncl
             % The final azimuth angles from each cluster i    
             phi_AOD = randLaplacian(Nray, 1, phiT(j,i), sig^2);
             phi_AOA = randLaplacian(Nray, 1, phiR(j,i), sig^2);    
             for l=1:Nray
                 at = exp(-1j*pi*sin(phi_AOD(l))*(0:1:NT-1)).'/sqrt(NT);
                 ar = exp(-1j*pi*sin(phi_AOA(l))*(0:1:NR-1)).'/sqrt(NR);
                 Alpha = alpha(j,l,i);
                 H = H + Alpha*ar*at'; % channel
                 At = [At at];
             end
         end
         H_n = sqrt(NT*NR/(Ncl*Nray))*H;
         ff(i) = norm(H_n,'fro')^2;
         [U,~,V] = svd(H_n);
         Fopt = V(:,1:NS);
%          Capacity_fully_digital(i) = log2(det(eye(NR)+(EN(k)/NS)*H_n*Fopt*Fopt'*H_n'));
         
         [FRF_mod,FBB_mod] = Mod_OMP(NS,NT_RF(m),NT/NT_RF(m),Fopt,At);
         Capacity_mod(i,k,m) = log2(det(eye(NR)+(EN(k)/NS)*H_n*FRF_mod*FBB_mod*FBB_mod'*FRF_mod'*H_n'));
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cap_mod(k,m) = mean(abs(Capacity_mod(:,k,m)));
    end
% cap_full_digital(k) = mean(abs(Capacity_fully_digital));
end

plot(SNR,smooth(cap_mod(:,1)),'-ks','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor',[0.1 0.1 0.1]);
hold on;
plot(SNR,smooth(cap_mod(:,2)),'-ys','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor','yellow','MarkerFaceColor',[1 .8 .0]);
hold on;
plot(SNR,smooth(cap_mod(:,3)),'-rs','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
hold on;
plot(SNR,smooth(cap_mod(:,4)),'-bs','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor','blue','MarkerFaceColor',[.3 .5 .8])
hold on;
plot(SNR,smooth(cap_mod(:,5)),'-gs','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor','green','MarkerFaceColor',[0 .5 .5]);
hold on;
plot(SNR,smooth(cap_mod(:,6)),'-cs','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor',	[0, 0.75, 0.75],'MarkerFaceColor',[0.3010, 0.7450, 0.9330]);
hold on;
plot(SNR,smooth(cap_mod(:,7)),'-ms','LineWidth',1,'MarkerSize',7,'MarkerEdgeColor',	[0.4940, 0.1840, 0.5560],'MarkerFaceColor',	[0.75, 0, 0.75]);
% hold on;
% plot(SNR,smooth(cap_full_digital),'-s','LineColor',	[0.8500, 0.3250, 0.0980],'LineWidth',1,'MarkerSize',7,'MarkerEdgeColor',	[0.4940, 0.1840, 0.5560],'MarkerFaceColor',	[0.6, 0.3, 0.6]);

xlim([minSNR maxSNR]);
ylim([minCap maxCap]);
grid on;
% title("Spectral efficiencies achieved by the proposed hybrid analog-digital precoding algorithm in a 64×16 mmWave system");

yticks(minCap:CapStep:maxCap);
set(get(gca,'XLabel'),'String','SNR(dB)','Interpreter','latex');
set(get(gca,'YLabel'),'String','Capacity (bps/Hz)','Interpreter','latex');
grid on;

h1 = legend('RF chain = 1',...
    'RF chain = 2',...
    'RF chain = 4',...
    'RF chain = 8',...
    'RF chain = 16',...
    'RF chain = 32',...
    'RF chain = 64',...
           'Location','Northwest');
set(h1, 'Fontsize', 10,'Interpreter','latex');