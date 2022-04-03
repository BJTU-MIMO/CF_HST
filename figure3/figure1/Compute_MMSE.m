function [v,u,A] = Compute_MMSE(h,hhat,hmean,C,I,ID,K,L,M,N,p,sigma2,nbrOfRealizations)

a = zeros(L,K,M);
test = zeros(N,N,nbrOfRealizations,K,L,K,M,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for l = 1:L
            for s = 1:M
                for i = 1:K
                    for m = 1:M
                        test(:,:,n,k,l,i,s,m) = (I(i,l,s,m)*hmean(:,n,i,l)+ID(s,m)*hhat(:,n,i,l))*(I(i,l,s,m)*hmean(:,n,i,l)+ID(s,m)*hhat(:,n,i,l))'+ID(s,m)^2*C(:,:,i,l);
                    end
                end
            end
        end
    end
end

test1 = squeeze(sum(p*sum(test(:,:,:,:,:,:,:,:),8),6));
v = zeros(N,nbrOfRealizations,K,L,M);
v1 = zeros(nbrOfRealizations,K,L,M);
u = zeros(L,nbrOfRealizations,K,K,M,M);

for n = 1:nbrOfRealizations
    for k = 1:K
        for l = 1:L
            for s = 1:M
                v(:,n,k,l,s) = p*((test1(:,:,n,k,l,s)+sigma2*eye(N))^(-1))*(I(k,l,s,s)*hmean(:,n,k,l)+ID(s,s)*hhat(:,n,k,l));
                v1(n,k,l,s) = norm(v(:,n,k,l,s))^2;
                for i = 1:K
                    for m = 1:M
                        u(l,n,k,i,s,m) = v(:,n,k,l,s)'*(I(i,l,s,m)*hmean(:,n,i,l)+ID(s,m)*h(:,n,i,l));
                    end
                end
            end
        end
    end
end

v2 = squeeze(mean(v1(:,:,:,:),1));
A = zeros(L,L,K,M);

for k = 1:K
    for l =1:L
        for s = 1:M
            A(l,l,k,s) = v2(k,l,s);
        end
    end
end


