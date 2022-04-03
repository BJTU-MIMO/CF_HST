function [SE] = Compute_SE_Cellular(ghat_c,hhat_c,hmean_c,C_c,I_c,ID_c,K,L,M,N,p,sigma2,nbrOfRealizations)

%a_kls = 1/L*ones(K,L,M);
%a_kls = repmat(beta_kl(:,:)./repmat(sum(beta_kl(:,:),1),[K 1]),[1 1 M]);

test = zeros(L*N,L*N,K,M,M);
term1 = zeros(nbrOfRealizations,K,M);
term2 = zeros(nbrOfRealizations,K,M,M);
term3 = zeros(nbrOfRealizations,K,K,M,M);
term4 = zeros(nbrOfRealizations,K,M);


for i = 1:K
    for s = 1:M
        for m = 1:M
            test(:,:,i,s,m) = ID_c(s,m)^2*C_c(:,:,i);
        end
    end
end

for n = 1:nbrOfRealizations
    for k = 1:K
        for s = 1:M
            term1(n,k,s) = ghat_c(:,n,k)'*(I_c(k,s,s)*hmean_c(:,n,k)+ID_c(s,s)*hhat_c(:,n,k));
            for m = 1:M
                if m == s
                    term2(n,k,s,m) = 0;
                else
                    term2(n,k,s,m) = ghat_c(:,n,k)'*(I_c(k,s,m)*hmean_c(:,n,k)+ID_c(s,m)*hhat_c(:,n,k));
                end
                for i = 1:K
                    if i == k
                        term3(n,k,i,s,m) = 0;
                    else
                        term3(n,k,i,s,m) = ghat_c(:,n,k)'*(I_c(i,s,m)*hmean_c(:,n,i)+ID_c(s,m)*hhat_c(:,n,i));
                    end
                end
            end
            term4(n,k,s) = ghat_c(:,n,k)'*(squeeze(sum(p*sum(test(:,:,:,s,:),5),3))+sigma2*eye(L*N))*ghat_c(:,n,k);
        end
    end
end

DS = p*abs(term1(:,:,:)).^2;
IT1 = squeeze(p*sum(abs(term2(:,:,:,:)).^2,4));
IT2 = squeeze(sum(p*sum(abs(term3(:,:,:,:,:)).^2,5),3));
IT3 = term4(:,:,:);

SINR = DS(:,:,:)./(IT1(:,:,:)+IT2(:,:,:)+IT3(:,:,:));
SE = squeeze(mean(log2(1+SINR(:,:,:)),1));



