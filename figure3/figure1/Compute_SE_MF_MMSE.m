function [SE] = Compute_SE_MF_MMSE(u,A,K,L,M,p,sigma2,nbrOfRealizations)

a= 1/L*ones(L,K,M);

test1 = zeros(K,M);
test2 = zeros(nbrOfRealizations,K,K,M,M);
test3 = zeros(K,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            test1(k,s) = a(:,k,s)'*squeeze(mean(u(:,:,k,k,s,s),2));
            for i = 1:K
                for m = 1:M
                    test2(n,k,i,s,m) = a(:,k,s)'*u(:,n,k,i,s,m);
                end
            end
            test3(k,s) = a(:,k,s)'*A(:,:,k,s)*a(:,k,s);
        end
    end
end

term1 = p*abs(test1(:,:)).^2;
term2 = squeeze(sum(p*sum(mean(abs(test2(:,:,:,:,:)).^2,1),5),3));
term3 = sigma2*test3(:,:);

SINR = term1(:,:)./(term2(:,:)-term1(:,:)+term3(:,:));

SE = log2(1+SINR(:,:));




