function [SE] = Compute_SE_SmallCell_MMSE(hhat,hmean,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations)


test1 = zeros(N,N,nbrOfRealizations,K,L,M,M);
test2 = zeros(N,N,nbrOfRealizations,K,K,L,M,M);
test3 = zeros(N,N,K,L,M,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for l = 1:L         
            for s = 1:M
                for m = 1:M
                    if m == s
                        test1(:,:,n,k,l,s,m) = 0;
                    else
                        test1(:,:,n,k,l,s,m) = (I(k,l,s,m)*hmean(:,n,k,l) + ID(s,m)*hhat(:,n,k,l))*(I(k,l,s,m)*hmean(:,n,k,l) + ID(s,m)*hhat(:,n,k,l))';
                    end
                    for i = 1:K
                        if i == k
                            test2(:,:,n,k,i,l,s,m) = 0;
                        else
                            test2(:,:,n,k,i,l,s,m) = (I(i,l,s,m)*hmean(:,n,i,l) + ID(s,m)*hhat(:,n,i,l))*(I(i,l,s,m)*hmean(:,n,i,l) + ID(s,m)*hhat(:,n,i,l))';
                        end
                        test3(:,:,i,l,s,m) = ID(s,m)^2*C(:,:,i,l);
                    end
                end
            end
        end
    end
end

term1 = squeeze(p*sum(test1(:,:,:,:,:,:,:),7));
term2 = squeeze(sum(p*sum(test2(:,:,:,:,:,:,:,:),8),5));
term3 = squeeze(sum(p*sum(test3(:,:,:,:,:,:),6),3));

SINR = zeros(nbrOfRealizations,K,L,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for l = 1:L
            for s = 1:M
                SINR(n,k,l,s) = p*(I(k,l,s,s)*hmean(:,n,k,l) + ID(s,s)*hhat(:,n,k,l))'*((term1(:,:,n,k,l,s) + term2(:,:,n,k,l,s) + term3(:,:,l,s) + sigma2*eye(N))^(-1))*(I(k,l,s,s)*hmean(:,n,k,l) + ID(s,s)*hhat(:,n,k,l));
            end
        end
    end
end

SE_kls = squeeze(mean(log2(1+SINR(:,:,:,:)),1));
SE = squeeze(max(SE_kls,[],2));


