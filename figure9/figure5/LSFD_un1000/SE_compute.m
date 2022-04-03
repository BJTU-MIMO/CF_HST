
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

d_ve = 50;

AP_place = linspace(-length_AP/2,length_AP/2,L);
UE_place = linspace(-length_train/2,length_train/2,K);

%distance = length_AP-length_train;
distance = 600;
Number = 600;
mobile = linspace(-distance/2,distance/2,Number);

speed = 100:50:600;    % km/h

ASDdeg = 10;

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

SE_n_LSFD = zeros(K,M,Number,length(speed));

for n = 1:Number
    disp(['n------------------- ' num2str(n) ' out of ' num2str(Number)]);
    for sp = 1:length(speed)
        disp(['s----------- ' num2str(sp) ' out of ' num2str(length(speed))]);
        
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
        eta_kil = zeros(K,K,L);
        
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
                RicianK = 1000;
                beta_LoS(k,l) = (RicianK/(RicianK+1))*beta(k,l);
                beta_NLoS(k,l) = (1/(RicianK+1))*beta(k,l);
                gamma(k,l) = p*tau_p*beta_NLoS(k,l)^2/(p*tau_p*beta_NLoS(k,l)+sigma2);
                
                hmean(:,k,l) = (sqrt(beta_LoS(k,l))*exp(1i*2*pi.*(0:(N-1))*sin_kl(k,l)*antennaSpacing)) * exp(1i*phase(k,l));
                correlationMatrix = functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);
                R(:,:,k,l) = beta_NLoS(k,l)*eyeN;
                %R(:,:,k,l) = beta_NLoS(k,l)*correlationMatrix;
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
        
        [SE_LSFD] = Compute_SE_LSFD(R,Q,hmean,I,ID,K,L,M,p,sigma2);
        SE_n_LSFD(:,:,n,sp) = SE_LSFD(:,:);
        
    end
end

%%

LSFD = squeeze(mean(mean(SE_n_LSFD(:,:,:,:),2),1));

%%

hold on; box on;

s = 1;
%Number_new = Number-2*s;
ss = (s+1):(Number-s);

plot(speed,squeeze(mean(LSFD(ss,:),1)),'r-o','LineWidth',2);

xlabel('Speed (km/h)','Interpreter','latex');
ylabel('Average SE (bit/s/Hz)','Interpreter','latex');
% legend('$K = 8$, $N = 4$, $d_{\mathrm{ve}} = 50$m','$K = 6$, $N = 2$, $d_{\mathrm{ve}} = 50$m','$K = 8$, $N = 2$, $d_{\mathrm{ve}} = 50$m','$K = 8$, $N = 2$, $d_{\mathrm{ve}} = 200$m','Interpreter','latex');
% set(gca, 'Fontname', 'Times New Roman','FontSize',14);

grid on;
%axis([100 600 0 5]);

