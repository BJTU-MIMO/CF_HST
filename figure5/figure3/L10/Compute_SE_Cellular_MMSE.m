function [SE] = Compute_SE_Cellular_MMSE(hhat_c,hmean_c,C_c,I_c,ID_c,K,L,M,N,p,sigma2,nbrOfRealizations)

test1 = zeros(L*N,L*N,nbrOfRealizations,K,M,M);
test2 = zeros(L*N,L*N,nbrOfRealizations,K,K,M,M);
test3 = zeros(L*N,L*N,K,M,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            for m = 1:M
                if m == s
                    test1(:,:,n,k,s,m) = 0;
                else
                    test1(:,:,n,k,s,m) = (I_c(k,s,m)*hmean_c(:,n,k) + ID_c(s,m)*hhat_c(:,n,k))*(I_c(k,s,m)*hmean_c(:,n,k) + ID_c(s,m)*hhat_c(:,n,k))';
                end
                for i = 1:K
                    if i == k
                        test2(:,:,n,k,i,s,m) = 0;
                    else
                        test2(:,:,n,k,i,s,m) = (I_c(i,s,m)*hmean_c(:,n,i) + ID_c(s,m)*hhat_c(:,n,i))*(I_c(i,s,m)*hmean_c(:,n,i) + ID_c(s,m)*hhat_c(:,n,i))';
                    end
                    test3(:,:,i,s,m) = ID_c(s,m)^2*C_c(:,:,i);
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
            SINR(n,k,s) = p*(I_c(k,s,s)*hmean_c(:,n,k) + ID_c(s,s)*hhat_c(:,n,k))'*((term1(:,:,n,k,s) + term2(:,:,n,k,s) + term3(:,:,s) + sigma2*eye(L*N))^(-1))*(I_c(k,s,s)*hmean_c(:,n,k) + ID_c(s,s)*hhat_c(:,n,k));
        end
    end
end
SE = squeeze(mean(log2(1+SINR(:,:,:)),1));

