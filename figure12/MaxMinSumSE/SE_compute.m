
clear;

%%

K = 8;
L = 20;
N = 4;
M = 2;

antennaSpacing = 1/2;

length_AP = 1000;
length_train = 200;
d_AP = length_AP/(L-1);
d_UE = length_train/(K-1);

d_ve = 10;

AP_place = linspace(-length_AP/2,length_AP/2,L);
UE_place = linspace(-length_train/2,length_train/2,K);

%distance = length_AP-length_train;
distance = 600;
Number = 600;
mobile = linspace(-distance/2,distance/2,Number);

speed = 300;    % km/h

ASDdeg = 30;

B = 20e6;
noiseFigure = 9;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
noiseVariance = db2pow(noiseVariancedBm)/1000; % W
sigma2 = noiseVariance;

p = 0.2;
tau_p = K;
eyeN = eye(N);

%threshold = -60:20:0;
threshold = -10;
beta_0 = -120;
ref = 1000;
aa = 3;

%%

% SE_n_MF_Scalable = zeros(K,M,Number,length(speed));
% SE_n_MF_Scalable_Fractional = zeros(K,M,Number,length(speed));
% SE_n_MF_Scalable_MaxMin = zeros(K,M,Number,length(speed));
% SE_n_MF_Scalable_SumSE = zeros(K,M,Number,length(speed));

SE_n_LSFD_Scalable = zeros(K,M,Number,length(speed));
SE_n_LSFD_Scalable_Fractional = zeros(K,M,Number,length(speed));
SE_n_MF_Scalable_MaxMin_LSFD = zeros(K,M,Number,length(speed));
SE_n_MF_Scalable_SumSE_LSFD = zeros(K,M,Number,length(speed));

for n = 1:Number %round(2*NumberOfSymbol/3);
    disp(['n------------------- ' num2str(n) ' out of ' num2str(Number)]);
    for sp = 1:length(speed)
        disp(['sp---- ' num2str(sp) ' out of ' num2str(length(speed))]);
        
        v = 1000*speed(sp)/3600;   % m/s
        T = 5 * 10^(-4);
        f = 2*10^9;
        c = 3*10^8;
        w = f*v*T/c;
        
        dho = zeros(K,L);
        d = zeros(K,L);
        cos_kl = zeros(K,L);
        sin_kl = zeros(K,L);
        epsilon = zeros(K,L);
        beta = zeros(K,L);
        beta_LoS = zeros(K,L);
        beta_NLoS = zeros(K,L);
        gamma = zeros(K,L);
        hmean = zeros(N,K,L);
        R = zeros(N,N,K,L);
        Psi = zeros(N,N,K,L);
        Q = zeros(N,N,K,L);
        C = zeros(N,N,K,L);
        I = zeros(K,L,M,M);
        ID = zeros(M,M);
        
        phase = - pi + (2*pi)*rand(K,L);
        
        for k = 1:K
            for l = 1:L
                dho(k,l) = (AP_place(l) - (mobile(n) + UE_place(k)));
                d(k,l) = (sqrt(dho(k,l)^2+d_ve^2));
                cos_kl(k,l) = (dho(k,l)/d(k,l));
                %angletoUE = acos(cos_kl(k,l));
                sin_kl(k,l) = cos_kl(k,l);
                angletoUE = asin(sin_kl(k,l));
                epsilon(k,l) = w*cos_kl(k,l);
                %beta(k,l) = (d(k,l)^(-4));
                beta(k,l) = db2pow(beta_0 + 10*log10((d(k,l)/ref)^(-aa)));
                %beta(k,l) = db2pow(-30.18 - 26*log10(d(k,l)));
                %beta(k,l) = db2pow(-50 - 10*log10(d(k,l)^2));
                RicianK = 100;
                beta_LoS(k,l) = (RicianK/(RicianK+1))*beta(k,l);
                beta_NLoS(k,l) = (1/(RicianK+1))*beta(k,l);
                gamma(k,l) = p*tau_p*beta_NLoS(k,l)^2/(p*tau_p*beta_NLoS(k,l)+sigma2);
                
                hmean(:,k,l) = (sqrt(beta_LoS(k,l))*exp(1i*2*pi.*(0:(N-1))*sin_kl(k,l)*antennaSpacing)) * exp(1i*phase(k,l));
                correlationMatrix = functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
                %R(:,:,k,l) = beta_NLoS(k,l)*eyeN;
                R(:,:,k,l) = beta_NLoS(k,l)*correlationMatrix;
                Psi(:,:,k,l) = (p*tau_p*R(:,:,k,l)+sigma2*eyeN)^(-1);
                Q(:,:,k,l) = p*tau_p*R(:,:,k,l)*Psi(:,:,k,l)*R(:,:,k,l);
                C(:,:,k,l) = R(:,:,k,l) - Q(:,:,k,l);
                
                for s = 1:M
                    for m = 1:M
                        I(k,l,s,m) = sin(pi*(m+epsilon(k,l)-s))/(M*sin(pi/M*(m+epsilon(k,l)-s)))*exp(1i*pi*(1-1/M)*(m+epsilon(k,l)-s));
                        if s == m
                            ID(s,m) = 1;
                        else
                            ID(s,m) = ((-1)^(m-s))*w/(sqrt(2)*(m-s));
                        end
                    end
                end
            end
            
        end
        
        D = zeros(K,L);
        %D = ones(K,L);
        pilotIndex = zeros(K,1);
        masterAPs = zeros(K,1);
        
        for k = 1:K
            [~,master] = max(beta(k,:));
            D(k,master) = 1;
            masterAPs(k) = master;
            if k <= tau_p
                pilotIndex(k) = k;
            else
                pilotinterference = zeros(tau_p,1);
                for t = 1:tau_p
                    pilotinterference(t) = sum(db2pow(beta(master,pilotIndex(1:k-1,n)==t)));
                end
                [~,bestpilot] = min(pilotinterference);
                pilotIndex(k) = bestpilot;
            end
        end
        
        for l = 1:L
            for t = 1:tau_p
                pilotUEs = find(t==pilotIndex(:));
                if sum(D(pilotUEs,l)) == 0
                    [gainValue,UEindex] = max(beta(pilotUEs,l));
                    if 10*log10(gainValue) - 10*log10(beta(pilotUEs(UEindex),masterAPs(pilotUEs(UEindex)))) >= threshold
                        D(pilotUEs(UEindex),l) = 1;
                    end
                end
            end
        end
        
        % min
        tt_min = 1;
        numl_min = squeeze(sum(D(:,:),2));
        zeta_i_min = (squeeze(sum(1./abs(d(:,:)).*D(:,:),2))./numl_min(:)).^(tt_min);
        p_min = min(zeta_i_min);
        eta_min = p_min./zeta_i_min(:);
        % max
        tt_max = 1;
        numl_max = squeeze(sum(D(:,:),2));
        zeta_i_max = (squeeze(sum(beta(:,:).*D(:,:),2))./numl_max(:)).^(tt_max);
        p_max = max(zeta_i_max);
        eta_max = zeta_i_max(:)/p_max;
        
        % MF -------------------------------------------------------------------------------------------
        %         [SE_MF_Scalable] = Compute_SE_MF_Scalable(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        %         SE_n_MF_Scalable(:,:,n,sp) = SE_MF_Scalable(:,:);
        %
        %         [SE_MF_Scalable_Fractional] = Compute_SE_MF_Scalable_Fractional(R,Q,hmean,I,ID,K,L,M,p,eta_min,sigma2,D);
        %         SE_n_MF_Scalable_Fractional(:,:,n,sp) = SE_MF_Scalable_Fractional(:,:);
        %
        %         [SE_MF_Scalable_MaxMin] = Compute_SE_MF_Scalable_MaxMin(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        %         SE_n_MF_Scalable_MaxMin(:,:,n,sp) = SE_MF_Scalable_MaxMin(:,:);
        %
        %         [SE_MF_Scalable_SumSE] = Compute_SE_MF_Scalable_SumSE(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        %         SE_n_MF_Scalable_SumSE(:,:,n,sp) = SE_MF_Scalable_SumSE(:,:);
        
        % LSFD -------------------------------------------------------------------------------------------
        
        [SE_LSFD_Scalable] = Compute_SE_LSFD_Scalable(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        SE_n_LSFD_Scalable(:,:,n,sp) = SE_LSFD_Scalable(:,:);
        
        [SE_LSFD_Scalable_Fractional] = Compute_SE_LSFD_Scalable_Fractional(R,Q,hmean,I,ID,K,L,M,p,eta_min,sigma2,D);
        SE_n_LSFD_Scalable_Fractional(:,:,n,sp) = SE_LSFD_Scalable_Fractional(:,:);
        
        [SE_MF_Scalable_MaxMin_LSFD] = Compute_SE_MF_Scalable_MaxMin_LSFD(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        SE_n_MF_Scalable_MaxMin_LSFD(:,:,n,sp) = SE_MF_Scalable_MaxMin_LSFD(:,:);
        
        [SE_MF_Scalable_SumSE_LSFD] = Compute_SE_MF_Scalable_SumSE_LSFD(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);
        SE_n_MF_Scalable_SumSE_LSFD(:,:,n,sp) = SE_MF_Scalable_SumSE_LSFD(:,:);
        
    end
end

%%

% MF_Scalable = squeeze(mean(mean(SE_n_MF_Scalable(:,:,:,:),2),1));
% MF_Scalable_Fractional = squeeze(mean(mean(SE_n_MF_Scalable_Fractional(:,:,:,:),2),1));
% MF_Scalable_MaxMin = squeeze(mean(mean(SE_n_MF_Scalable_MaxMin(:,:,:,:),2),1));
% MF_Scalable_SumSE = squeeze(mean(mean(SE_n_MF_Scalable_SumSE(:,:,:,:),2),1));

% LSFD_Scalable = squeeze(mean(mean(SE_n_LSFD_Scalable(:,:,:,:),2),1));
% LSFD_Scalable_Fractional = squeeze(mean(mean(SE_n_LSFD_Scalable_Fractional(:,:,:,:),2),1));
% MF_Scalable_MaxMin_LSFD = squeeze(mean(mean(SE_n_MF_Scalable_MaxMin_LSFD(:,:,:,:),2),1));
% MF_Scalable_SumSE_LSFD = squeeze(mean(mean(SE_n_MF_Scalable_SumSE_LSFD(:,:,:,:),2),1));

%%
subplot(2,1,1);
hold on; box on;

s = 0;
Number_new = Number-2*s;
ss = (s+1):(Number-s);

t = 1;

% plot(sort(reshape(SE_n_MF_Scalable(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'r-','LineWidth',2);
% plot(sort(reshape(SE_n_MF_Scalable_Fractional(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'b--','LineWidth',2);
% plot(sort(reshape(SE_n_MF_Scalable_MaxMin(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'k-.','LineWidth',2);
% plot(sort(reshape(SE_n_MF_Scalable_SumSE(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'m:','LineWidth',2);

plot(sort(reshape(SE_n_LSFD_Scalable(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'r-','LineWidth',2);
plot(sort(reshape(SE_n_LSFD_Scalable_Fractional(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'b--','LineWidth',2);
plot(sort(reshape(SE_n_MF_Scalable_MaxMin_LSFD(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'k-.','LineWidth',2);
plot(sort(reshape(SE_n_MF_Scalable_SumSE_LSFD(:,:,:,t),[K*M*Number 1])),linspace(0,1,K*M*Number),'m:','LineWidth',2);

% plot(speed,squeeze(mean(MF_Scalable(ss,:),1)),'r-o','LineWidth',2);
% plot(speed,squeeze(mean(MF_Scalable_Fractional(ss,:),1)),'b--<','LineWidth',2);
% plot(speed,squeeze(mean(MF_Scalable_MaxMin(ss,:),1)),'k-.+','LineWidth',2);
% plot(speed,squeeze(mean(MF_Scalable_SumSE(ss,:),1)),'m:*','LineWidth',2);

% plot(speed,squeeze(mean(LSFD_Scalable(ss,:),1)),'r-o','LineWidth',2);
% plot(speed,squeeze(mean(LSFD_Scalable_Fractional(ss,:),1)),'b--<','LineWidth',2);
% plot(speed,squeeze(mean(MF_Scalable_MaxMin_LSFD(ss,:),1)),'k-.+','LineWidth',2);
% plot(speed,squeeze(mean(MF_Scalable_SumSE_LSFD(ss,:),1)),'m:*','LineWidth',2);

grid on;
%axis([2 6 0 1]);


