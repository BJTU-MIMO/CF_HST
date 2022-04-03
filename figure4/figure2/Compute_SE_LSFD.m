function [SE] = Compute_SE_LSFD(R,Q,hmean,I,ID,K,L,M,p,sigma2)

a_kls = 1/L*ones(K,L,M);
%a_kls = repmat(beta_kl(:,:)./repmat(sum(beta_kl(:,:),1),[K 1]),[1 1 M]);

a = zeros(L,K,M);
b = zeros(L,K);
T = zeros(L,L,K,K,M,M);
c = zeros(L,K,M,M);
d = zeros(L,K,K,M,M);
A = zeros(L,L,K);

for l = 1:L
    for k = 1:K
        b(l,k) = I(k,l,1,1)*hmean(:,k,l)'*hmean(:,k,l) + ID(1,1)*trace(Q(:,:,k,l));
        for i = 1:K
            for s = 1:M
                a(l,k,s) = a_kls(k,l,s);
                for m = 1:M
                    T(l,l,k,i,s,m) = ID(s,m)^2*trace(R(:,:,i,l)*Q(:,:,k,l)) + (ID(s,m)^2*hmean(:,k,l)'*R(:,:,i,l)*hmean(:,k,l)) + I(i,l,s,m)^2*hmean(:,i,l)'*Q(:,:,k,l)*hmean(:,i,l);
                    c(l,k,s,m) = I(k,l,s,m)*hmean(:,k,l)'*hmean(:,k,l) + ID(s,m)*trace(Q(:,:,k,l));
                    d(l,k,i,s,m) = I(i,l,s,m)*hmean(:,k,l)'*hmean(:,i,l);                    
                end
            end
        end
        A(l,l,k) = trace(Q(:,:,k,l)) + hmean(:,k,l)'*hmean(:,k,l);
    end
end

term1 = zeros(L,L,K,M,M);
term2 = zeros(L,L,K,K,M,M);
SINR = zeros(K,M);

for k = 1:K
    for s = 1:M
        for i = 1:K
            for m = 1:M
                if m == s
                    term1(:,:,k,s,m) = 0;
                else
                    term1(:,:,k,s,m) = c(:,k,s,m)*c(:,k,s,m)';
                end
                if i == k
                    term2(:,:,k,i,s,m) = 0;
                else
                    term2(:,:,k,i,s,m) = d(:,k,i,s,m)*d(:,k,i,s,m)';
                end
            end
        end
        SINR(k,s) = p*b(:,k)'*(squeeze(sum(p*sum(T(:,:,k,:,s,:),6),4))+squeeze(p*sum(term1(:,:,k,s,:),5))+squeeze(sum(p*sum(term2(:,:,k,:,s,:),6),4))+sigma2*A(:,:,k))^(-1)*b(:,k);
    end
end

SE = log2(1+SINR(:,:));

