function [SE] = Compute_SE_SmallCell(ghat,hhat,hmean,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations)

%a_kls = 1/L*ones(K,L,M);
%a_kls = repmat(beta_kl(:,:)./repmat(sum(beta_kl(:,:),1),[K 1]),[1 1 M]);

test = zeros(N,N,K,L,M,M);
term1 = zeros(nbrOfRealizations,K,L,M);
term2 = zeros(nbrOfRealizations,K,L,M,M);
term3 = zeros(nbrOfRealizations,K,K,L,M,M);
term4 = zeros(nbrOfRealizations,K,L,M);


for i = 1:K
    for l = 1:L
        for s = 1:M
            for m = 1:M
                test(:,:,i,l,s,m) = ID(s,m)^2*C(:,:,i,l);
            end
        end
    end
end

for n = 1:nbrOfRealizations
    for k = 1:K
        for l = 1:L         
            for s = 1:M
                 term1(n,k,l,s) = ghat(:,n,k,l)'*(I(k,l,s,s)*hmean(:,n,k,l)+ID(s,s)*hhat(:,n,k,l));
                for m = 1:M
                    if m == s
                       term2(n,k,l,s,m) = 0;
                    else
                       term2(n,k,l,s,m) = ghat(:,n,k,l)'*(I(k,l,s,m)*hmean(:,n,k,l)+ID(s,m)*hhat(:,n,k,l));
                    end
                    for i = 1:K
                        if i == k
                           term3(n,k,i,l,s,m) = 0;
                        else
                           term3(n,k,i,l,s,m) = ghat(:,n,k,l)'*(I(i,l,s,m)*hmean(:,n,i,l)+ID(s,m)*hhat(:,n,i,l));
                        end
                    end
                end
                term4(n,k,l,s) = ghat(:,n,k,l)'*(squeeze(sum(p*sum(test(:,:,:,l,s,:),6),3))+sigma2*eye(N))*ghat(:,n,k,l);
            end            
        end
    end
end

DS = p*abs(term1(:,:,:,:)).^2;
IT1 = squeeze(p*sum(abs(term2(:,:,:,:,:)).^2,5));
IT2 = squeeze(sum(p*sum(abs(term3(:,:,:,:,:,:)).^2,6),3));
IT3 = term4(:,:,:,:);

SINR = DS(:,:,:,:)./(IT1(:,:,:,:)+IT2(:,:,:,:)+IT3(:,:,:,:));
SE_kls = squeeze(mean(log2(1+SINR(:,:,:,:)),1));
SE = squeeze(max(SE_kls,[],2));


