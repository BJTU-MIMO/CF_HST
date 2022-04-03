function [SE] = Compute_SE_LSFD_MMSE(u,A,K,L,M,p,sigma2,nbrOfRealizations)

test1 = zeros(L,L,nbrOfRealizations,K,K,M,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            for i = 1:K
                for m = 1:M
                    test1(:,:,n,k,i,s,m) = u(:,n,k,i,s,m)*u(:,n,k,i,s,m)';
                end
            end
        end
    end
end

test2 = squeeze(sum(p*sum(mean(test1(:,:,:,:,:,:,:),3),7),5));
SINR = zeros(K,M);

for k = 1:K
    for s = 1:M
        SINR(k,s) = p*squeeze(mean(u(:,:,k,k,s,s),2))'*((test2(:,:,k,s)+sigma2*A(:,:,k,s)-p*squeeze(mean(u(:,:,k,k,s,s),2))*squeeze(mean(u(:,:,k,k,s,s),2))')^(-1))*squeeze(mean(u(:,:,k,k,s,s),2));      
    end
end

SE = log2(1+SINR(:,:));

