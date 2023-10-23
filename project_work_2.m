close all;
warning('off','all');

% bitstream length
T = 10^5;
MC = 1000;

% model parameters
K = 4;
L = 62;
N = 3;
R = 4;
lnr = L*N*R;

% snr probed
snr_db = -40:2:20;
snr = 10.^(snr_db/10);
vari = 1./snr;

% variables for code
Pe_sim = zeros(MC,length(snr_db));
Pe_sim_als = zeros(MC,length(snr_db));
Pe_sim_rals = zeros(MC,length(snr_db));
Pe_sim_zf = zeros(MC,length(snr_db));

% processing at CPU
% Monte Carlo 
for sim_ind = 1:1:MC
    disp(sim_ind)

    % get full rank X_mod
    Xb = 2*(rand(1,K*T)>0.5)-1;
    X_mod = reshape(Xb,T,K);

    % modelling H
    H = location(L,N,K);

    % get full rank D
    HdSrc = hadamard(R);
    D = HdSrc(:,1:K);
    
    % get full rank khatris
    while (rank(X_mod) ~= K) && (rank(H) ~= K) && (rank(khatrirao(D,H))~=K) && (rank(khatrirao(X_mod,D))~=K) && (rank(khatrirao(H,X_mod))~=K)
        Xb = 2*(rand(1,K*T)>0.5)-1;
        X_mod = reshape(Xb,T,K);
        H = location(L,N,K);
    end
 
    % short data for ZF
    Xb_z = reshape(X_mod(L*N+1:T, :),K*(T-L*N),1);

    % calculate Gamma
    G1 = khatrirao(D,H).';
    G2 = khatrirao(D,X_mod).';

    for snr_ind = 1:length(snr_db)
        disp(snr_ind)
        % noise formation
        V = sqrt(1/(2 * snr(snr_ind)))*(randn(L*N,T,R) + 1i*randn(L*N,T,R));

        % system model
        % tensor formation
        Y = full(ktensor({H,X_mod,D})) + V;
        Y1 = tenmat(Y,1).data;
        Y2 = tenmat(Y,2).data;
        Y_eqv = Y2;

        % detection and estimation
        % Non-blind detection (H,D known)
        X_dem_eqv = Y_eqv*pinv(G1);
        
        % ALS detection (D known)
        X_dem_als = ALS(L*N,T,K,Y1,Y2,D,X_mod(1,:));

        % RALS detection (D known)
        X_dem_rals = R_ALS(L*N,T,K,Y1,Y2,D,X_mod(1,:));

        % pilot blind detection (D known)
        X_dem_zf = ZF(L*N,T,Y,D,X_mod,V);

        % reshape
        Xb_dem = reshape(X_dem_eqv,1,K*T);
        Xb_dem_als = reshape(X_dem_als,1,K*T);
        Xb_dem_rals = reshape(X_dem_rals,1,K*T);
        Xb_dem_zf = reshape(X_dem_zf,1,K*(T-L*N));
        %Xb_dem_cd = reshape(X_dem_cd,1,K*T);

        % quantization
        Xb_hat_eqv = 2*(real(Xb_dem) > 0)-1;
        Xb_hat_als = 2*(real(Xb_dem_als) > 0)-1;
        Xb_hat_rals = 2*(real(Xb_dem_rals) > 0)-1;
        Xb_hat_zf = 2*(real(Xb_dem_zf) > 0)-1;
        %Xb_dem_cd = 2*(real(Xb_dem_cd) > 0)-1;

        % error calculation
        error_eqv = (norm(Xb-Xb_hat_eqv)^2)/4;
        error_als = (norm(Xb-Xb_hat_als)^2)/4;
        error_rals = (norm(Xb-Xb_hat_rals)^2)/4;
        error_zf = (norm(Xb_z-Xb_hat_zf)^2)/4;
        %error_cd = (norm(Xb-Xb_hat_cd)^2)/4;

        % simulation probability
        Pe_sim(sim_ind,snr_ind) = error_eqv/(K*T);
        Pe_sim_als(sim_ind,snr_ind) = error_als/(K*T);
        Pe_sim_rals(sim_ind,snr_ind) = error_rals/(K*T);
        Pe_sim_zf(sim_ind,snr_ind) = error_zf/(K*(T-L*N));
        %Pe_sim_cd(sim_ind,snr_ind) = error_cd/(K*T);
    end
end
Avg_Pe_sim = sum(Pe_sim)/MC;
Avg_Pe_als = sum(Pe_sim_als)/MC;
Avg_Pe_rals = sum(Pe_sim_rals)/MC;
Avg_Pe_zf = sum(Pe_sim_zf)/MC;
%Avg_Pe_cd = sum(Pe_sim_cd)/MC;

% Theoretical
dof = lnr - K + 1;
PeEx   = zeros(1,length(snr_db));
for snr_ind = 1:length(snr_db)
    SNR = 10^(snr_db(snr_ind)/10);
    mu  = sqrt(SNR/(1+SNR));
    prod_term = ((1-mu)/2)^dof;
    sumterm = 0;
    for l = 0:dof
        sumterm = sumterm + nchoosek((dof-1+l),l)*(((1+mu)/2)^l);
    end
    PeEx(snr_ind) = prod_term*sumterm;
end

% plots
semilogy(snr_db,Avg_Pe_sim,"linewidth",1.5);
hold on;
semilogy(snr_db,Avg_Pe_als,"linewidth",1.5);
hold on;
semilogy(snr_db,Avg_Pe_rals,"linewidth",1.5);
hold on;
semilogy(snr_db,Avg_Pe_zf,"linewidth",1.5);
hold on;
semilogy(snr_db,PeEx,'--',"linewidth",1.5);

% limits
xlim([-40 20]);
ylim([10^-5 1]);
xlabel('SNR (dB)');
ylabel('Pe');
title('BER vs SNR (K=4)');
grid on;
legend('Non-blind ZF','ALS','RALS','ZF','Theo');
saveas(gcf,'final_plot.png');