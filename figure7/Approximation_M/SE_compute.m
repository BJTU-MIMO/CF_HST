
clear;

%%

K = 8;
L = 20;
N = 4;
M = 1024;

antennaSpacing = 1/2;

length_AP = 1000;
length_train = 200;
d_AP = length_AP/(L-1);
d_UE = length_train/(K-1);

%d_ve = 5:15:100;
%d_ve = 5:45:230;
d_ve = 50;

AP_place = linspace(-length_AP/2,length_AP/2,L);
UE_place = linspace(-length_train/2,length_train/2,K);

speed = 300;    % km/h

v = 1000*speed/3600;   % m/s
T = 5 * 10^(-4);
f = 2*10^9;
c = 3*10^8;
w = f*v*T/c;

ASDdeg = 30;

B = 20e6;
noiseFigure = 9;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
noiseVariance = db2pow(noiseVariancedBm)/1000; % W
sigma2 = noiseVariance;

p = 0.2;
tau_p = K;
eyeN = eye(N);
eyeLN = eye(L*N);

beta_0 = -120;
ref = 1000;
aa = 3;

%%

SE_n_LSFD = zeros(K,length(M));

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

eta_kil = zeros(K,K,L);

phase = - pi + (2*pi)*rand(K,L);

for k = 1:K
    for l = 1:L
        dho(k,l) = (AP_place(l) - UE_place(k));
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
        
    end   
end

for sp = 1:length(M)
    disp(['s----------- ' num2str(sp) ' out of ' num2str(length(M))]);
    I = zeros(K,L,M(sp),M(sp));
    ID = zeros(M(sp),M(sp));
    for k = 1:K
       disp(['k1-- ' num2str(k) ' out of ' num2str(K)]);
        for l = 1:L
            for s = 1:M(sp)
                for m = 1:M(sp)
                    I(k,l,s,m) = sin(pi*(m+epsilon(k,l)-s))/(M(sp)*sin(pi/M(sp)*(m+epsilon(k,l)-s)))*exp(1i*pi*(1-1/M(sp))*(m+epsilon(k,l)-s));
                    if s == m
                        ID(s,m) = 1;
                    else
                        ID(s,m) = ((-1)^(m-s))*w/(sqrt(2)*(m-s));
                    end
                end
            end
        end
    end
    
    [SE_LSFD] = Compute_SE_LSFD(R,Q,hmean,I,ID,K,L,M(sp),p,sigma2);
    SE_n_LSFD(:,sp) = squeeze(mean(SE_LSFD(:,:),2));
    
end


%%

LSFD = squeeze(mean(SE_n_LSFD(:,:),1));

%%

hold on; box on;

plot(M,LSFD,'r-o','LineWidth',2);

xlabel('Speed (km/h)','Interpreter','latex');
ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
% legend('$K = 8$, $N = 4$, $d_{\mathrm{ve}} = 50$m','$K = 6$, $N = 2$, $d_{\mathrm{ve}} = 50$m','$K = 8$, $N = 2$, $d_{\mathrm{ve}} = 50$m','$K = 8$, $N = 2$, $d_{\mathrm{ve}} = 200$m','Interpreter','latex');
% set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
%axis([0 1024 1.5 3]);

