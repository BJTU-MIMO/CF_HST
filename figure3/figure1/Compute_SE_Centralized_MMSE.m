function [SE] = Compute_SE_Centralized_MMSE(ghat,hhat,hmean,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations)

%a_kls = 1/L*ones(K,L,M);
%a_kls = repmat(beta_kl(:,:)./repmat(sum(beta_kl(:,:),1),[K 1]),[1 1 M]);

ghat_LN = zeros(L*N,nbrOfRealizations,K);
hhat_LN = zeros(L*N,nbrOfRealizations,K);
hmean_LN = zeros(L*N,nbrOfRealizations,K);
C_LN = zeros(L*N,L*N,K);
I_LN = zeros(L*N,L*N,K,M,M);

for l = 1:L
    for k = 1:K
        ghat_LN((l-1)*N+1:l*N,:,k) = ghat(:,:,k,l);
        hhat_LN((l-1)*N+1:l*N,:,k) = hhat(:,:,k,l);
        hmean_LN((l-1)*N+1:l*N,:,k) = hmean(:,:,k,l);
        C_LN((l-1)*N+1:l*N,(l-1)*N+1:l*N,k) = C(:,:,k,l);
        for s = 1:M
            for m = 1:M
                I_LN((l-1)*N+1:l*N,(l-1)*N+1:l*N,k,s,m) = I(k,l,s,m)*eye(N);
            end
        end
    end
end

test1 = zeros(L*N,L*N,nbrOfRealizations,K,M,M);
test2 = zeros(L*N,L*N,nbrOfRealizations,K,K,M,M);
test3 = zeros(L*N,L*N,nbrOfRealizations,K,M,M);


for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            for m = 1:M
                %f_LN(:,n,k,s,m) = I_LN(:,:,k,s,m)*hmean_LN(:,n,k)+ID(s,m)*hhat_LN(:,n,k);
                if m == s
                    test1(:,:,n,k,s,m) = 0;
                else
                    test1(:,:,n,k,s,m) = (I_LN(:,:,k,s,m)*hmean_LN(:,n,k)+ID(s,m)*hhat_LN(:,n,k))*(I_LN(:,:,k,s,m)*hmean_LN(:,n,k)+ID(s,m)*hhat_LN(:,n,k))'; 
                end
                for i = 1:K
                    if i == k
                       test2(:,:,n,k,i,s,m) = 0;
                    else
                       test2(:,:,n,k,i,s,m) = (I_LN(:,:,i,s,m)*hmean_LN(:,n,i)+ID(s,m)*hhat_LN(:,n,i))*(I_LN(:,:,i,s,m)*hmean_LN(:,n,i)+ID(s,m)*hhat_LN(:,n,i))';  
                    end
                       test3(:,:,i,s,m) = ID(s,m)^2*C_LN(:,:,i);
                end
            end
        end
    end
end

term1 = squeeze(p*sum(test1(:,:,:,:,:,:),6));
term2 = squeeze(sum(p*sum(test2(:,:,:,:,:,:,:),7),5));
term3 = squeeze(sum(p*sum(test3(:,:,:,:,:),5),3));
SINR = zeros(nbrOfRealizations,K,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            %SINR(n,k,s) = p*(I_LN(:,:,k,s,s)*hmean_LN(:,n,k)+ID(s,s)*hhat_LN(:,n,k))'*((term1(:,:,n,k,s)+term2(:,:,n,s)+term3(:,:,s)+sigma2*eye(L*N))^(-1))*(I_LN(:,:,k,s,s)*hmean_LN(:,n,k)+ID(s,s)*hhat_LN(:,n,k));
            SINR(n,k,s) = p*(I_LN(:,:,k,s,s)*hmean_LN(:,n,k)+ID(s,s)*hhat_LN(:,n,k))'*((term1(:,:,n,k,s)+term2(:,:,n,k,s)+term3(:,:,s)+sigma2*eye(L*N))^(-1))*(I_LN(:,:,k,s,s)*hmean_LN(:,n,k)+ID(s,s)*hhat_LN(:,n,k));
        end
    end
end
SE = squeeze(mean(log2(1+SINR(:,:,:)),1));





