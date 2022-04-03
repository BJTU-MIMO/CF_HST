
clear;

%%

K = 8;
L = 10;
N = 4;
M = 2;

nbrOfRealizations = 50;
antennaSpacing = 1/2;

length_AP = 1000;
length_train = 200;
d_AP = length_AP/(L-1);
d_UE = length_train/(K-1);

d_ve = 50;

AP_place = linspace(-length_AP/2,length_AP/2,L);
UE_place = linspace(-length_train/2,length_train/2,K);

%distance = length_AP-length_train;
distance = 600;
Number = 600;
mobile = linspace(-distance/2,distance/2,Number);

speed = 300;    % km/h
v = 1000*speed/3600;   % m/s
T = 5 * 10^(-4);

f = 2*10^9;
c = 3*10^8;
w = f*v*T/c;

ASDdeg = 30;

B = 20e6;
%B = 15e3;
noiseFigure = 9;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
noiseVariance = db2pow(noiseVariancedBm)/1000; % W
sigma2 = noiseVariance;
Np = sqrt(0.5*sigma2)*(randn(N,nbrOfRealizations,K,L) + 1i*randn(N,nbrOfRealizations,K,L));
LNp = sqrt(0.5*sigma2)*(randn(L*N,nbrOfRealizations,K) + 1i*randn(L*N,nbrOfRealizations,K));

p = 0.2;
tau_p = K;
eyeN = eye(N);
eyeLN = eye(L*N);

beta_0 = -120;
ref = 1000;
aa = 3;

%%

SE_n_Centralized_MMSE = zeros(K,M,Number);
%SE_n_MF_MMSE = zeros(K,M,Number);
SE_n_LSFD_MMSE = zeros(K,M,Number);

for n = 1:Number
    disp(['n------------------- ' num2str(n) ' out of ' num2str(Number)]);
    
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
    eta_kil = zeros(K,K,L);
    
    h = zeros(N,nbrOfRealizations,K,L);
    g = zeros(N,nbrOfRealizations,K,L);
    z = zeros(N,nbrOfRealizations,K,L);
    hhat = zeros(N,nbrOfRealizations,K,L);
    ghat = zeros(N,nbrOfRealizations,K,L);
    hmean1 = zeros(N,nbrOfRealizations,K,L);
    
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
            
            % channel estimation
            h(:,:,k,l) = sqrtm(R(:,:,k,l))*sqrt(0.5)*(randn(N,nbrOfRealizations)+1i*randn(N,nbrOfRealizations));
            %h(:,:,k,l) = sqrt(beta_NLoS(k,l))*sqrt(0.5)*(randn(N,nbrOfRealizations)+1i*randn(N,nbrOfRealizations));
            g(:,:,k,l) = h(:,:,k,l) + reshape(repmat(hmean(:,k,l),nbrOfRealizations,1),N,nbrOfRealizations);
            z(:,:,k,l) = sqrt(p*tau_p)*h(:,:,k,l)+Np(:,:,k,l);
            hhat(:,:,k,l) = sqrt(p*tau_p)*R(:,:,k,l)*Psi(:,:,k,l)*z(:,:,k,l);
            ghat(:,:,k,l) = hhat(:,:,k,l) + reshape(repmat(hmean(:,k,l),nbrOfRealizations,1),N,nbrOfRealizations);
            hmean1(:,:,k,l) = reshape(repmat(hmean(:,k,l),nbrOfRealizations,1),N,nbrOfRealizations);
            
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
    
    [SE_Centralized_MMSE] = Compute_SE_Centralized_MMSE(ghat,hhat,hmean1,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations);
    SE_n_Centralized_MMSE(:,:,n) = SE_Centralized_MMSE(:,:);    
    
    %---------------------
    [v,u,A] = Compute_MMSE(h,hhat,hmean1,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations);
    
%     [SE_MF_MMSE] = Compute_SE_MF_MMSE(u,A,K,L,M,p,sigma2,nbrOfRealizations);
%     SE_n_MF_MMSE(:,:,n) = SE_MF_MMSE(:,:);
    
    [SE_LSFD_MMSE] = Compute_SE_LSFD_MMSE(u,A,K,L,M,p,sigma2,nbrOfRealizations);
    SE_n_LSFD_MMSE(:,:,n) = SE_LSFD_MMSE(:,:);
    
end

%%

Centralized_MMSE = squeeze(mean(mean(SE_n_Centralized_MMSE(:,:,:),2),1));
%MF_MMSE = squeeze(mean(mean(SE_n_MF_MMSE(:,:,:),2),1));
LSFD_MMSE = squeeze(mean(mean(SE_n_LSFD_MMSE(:,:,:),2),1));

%%
hold on; box on;

plot(mobile(1:Number),Centralized_MMSE(1:Number),'r-','LineWidth',1);
%plot(mobile(1:Number),MF_MMSE(1:Number),'b--','LineWidth',1);
plot(mobile(1:Number),LSFD_MMSE(1:Number),'k-.','LineWidth',1);

%axis([-300 300 0 11]);
%% CDF
hold on; box on;
plot(sort(reshape(SE_n_Centralized_MMSE(:,:,:),[K*M*Number 1])),linspace(0,1,K*M*Number),'r-','LineWidth',2);
%plot(sort(reshape(SE_n_MF_MMSE(:,:,:),[K*M*Number 1])),linspace(0,1,K*M*Number),'b--','LineWidth',2);
plot(sort(reshape(SE_n_LSFD_MMSE(:,:,:),[K*M*Number 1])),linspace(0,1,K*M*Number),'k-.','LineWidth',2);

%  axis([0 3 0 1]);



