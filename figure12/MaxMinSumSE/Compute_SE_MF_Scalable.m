 function [SE] = Compute_SE_MF_Scalable(R,Q,hmean,I,ID,K,L,M,p,sigma2,D)

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
        b(l,k) = I(k,l,1,1)*hmean(:,k,l)'*D(k,l)*hmean(:,k,l) + ID(1,1)*trace(D(k,l)*Q(:,:,k,l));
        for i = 1:K
            for s = 1:M
                a(l,k,s) = a_kls(k,l,s);
                for m = 1:M
                    T(l,l,k,i,s,m) = ID(s,m)^2*trace(R(:,:,i,l)*D(k,l)*Q(:,:,k,l)) + (ID(s,m)^2*hmean(:,k,l)'*D(k,l)*R(:,:,i,l)*D(k,l)*hmean(:,k,l)) + I(i,l,s,m)^2*hmean(:,i,l)'*D(k,l)*Q(:,:,k,l)*D(k,l)*hmean(:,i,l);
                    c(l,k,s,m) = I(k,l,s,m)*hmean(:,k,l)'*D(k,l)*hmean(:,k,l) + ID(s,m)*trace(D(k,l)*Q(:,:,k,l));
                    d(l,k,i,s,m) = I(i,l,s,m)*hmean(:,k,l)'*D(k,l)*hmean(:,i,l);                    
                end
            end
        end
        A(l,l,k) = trace(D(k,l)*Q(:,:,k,l)*D(k,l)) + hmean(:,k,l)'*D(k,l)*D(k,l)*hmean(:,k,l);
    end
end

term1 = zeros(K,M);
term2 = zeros(K,K,M,M);
term3 = zeros(K,M,M);
term4 = zeros(K,K,M,M);
term5 = zeros(K,M);

for k = 1:K
    for s = 1:M
        term1(k,s) = a(:,k,s)'*b(:,k);
        for i = 1:K
            for m = 1:M
                term2(k,i,s,m) = a(:,k,s)'*T(:,:,k,i,s,m)*a(:,k,s);
                if m == s
                   term3(k,s,m) = 0;
                else
                   term3(k,s,m) = a(:,k,s)'*c(:,k,s,m);
                end
                if i == k
                   term4(k,i,s,m) = 0;
                else
                   term4(k,i,s,m) = a(:,k,s)'*d(:,k,i,s,m);
                end
            end
        end
        term5(k,s) = a(:,k,s)'*A(:,:,k)*a(:,k,s);
    end
end

DS = p*abs(term1(:,:)).^2;
IT1 = squeeze(sum(p*sum(real(term2(:,:,:,:)),4),2));
IT2 = squeeze(p*sum(abs(term3(:,:,:)).^2,3));
IT3 = squeeze(sum(p*sum(abs(term4(:,:,:,:)).^2,4),2));
IT4 = sigma2*real(term5(:,:));

SINR = DS(:,:)./(IT1(:,:)+IT2(:,:)+IT3(:,:)+IT4(:,:));
SE = log2(1+SINR(:,:));

