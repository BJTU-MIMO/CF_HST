function [SE] = Compute_SE_LSFD_Scalable_Power(R,Q,hmean,I,ID,K,L,M,p,eta,sigma2,D)

% DD = round(D);
% [~,l_first] = find(DD(1,:),1,'first'); 
% [~,l_last] = find(DD(K,:),1,'last'); 
% 
% L_size = l_last - l_first + 1;

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
                for m = 1:M
                    T(l,l,k,i,s,m) = p*eta(i)*(ID(s,m)^2*trace(R(:,:,i,l)*D(k,l)*Q(:,:,k,l)) + (ID(s,m)^2*hmean(:,k,l)'*D(k,l)*R(:,:,i,l)*D(k,l)*hmean(:,k,l)) + I(i,l,s,m)^2*hmean(:,i,l)'*D(k,l)*Q(:,:,k,l)*D(k,l)*hmean(:,i,l));
                    c(l,k,s,m) = sqrt(p*eta(k))*(I(k,l,s,m)*hmean(:,k,l)'*D(k,l)*hmean(:,k,l) + ID(s,m)*trace(D(k,l)*Q(:,:,k,l)));
                    d(l,k,i,s,m) = sqrt(p*eta(i))*(I(i,l,s,m)*hmean(:,k,l)'*D(k,l)*hmean(:,i,l));                    
                end
            end
        end
        A(l,l,k) = trace(D(k,l)*Q(:,:,k,l)*D(k,l)) + hmean(:,k,l)'*D(k,l)*D(k,l)*hmean(:,k,l);
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
        
%         %b(all(b==0,2),k)=[];
%         G1 = squeeze(sum(p*sum(T(:,:,k,:,s,:),6),4));
%         %G1(find(all(G1==0,2)),:)=[];
%         %G1(:,find(all(G1==0,1))) = [];
%         G2 = squeeze(p*sum(term1(:,:,k,s,:),5));
%        % G2(find(all(G2==0,2)),:)=[];
%         %G2(:,find(all(G2==0,1))) = [];
%         G3 = squeeze(sum(p*sum(term2(:,:,k,:,s,:),6),4));
%         %G3(find(all(G3==0,2)),:)=[];
%         %G3(:,find(all(G3==0,1))) = [];
%         G4 = sigma2*A(:,:,k);
%         %G4(find(all(G4==0,2)),:)=[];
%         %G4(:,find(all(G4==0,1))) = [];
        
        %SINR(k,s) = p*b(:,k)'*(G1(:,:)+G2(:,:)+G4(:,:))^(-1)*b(:,k);
        bb = b(:,k);
        G = squeeze(sum(sum(T(:,:,k,:,s,:),6),4)) + squeeze(sum(term1(:,:,k,s,:),5)) + squeeze(sum(sum(term2(:,:,k,:,s,:),6),4)) + sigma2*A(:,:,k);
        %G = G3(:,:,k,s);
        bb(all(bb==0,2),:)=[];
        G(find(all(G==0,2)),:)=[];
        G(:,find(all(G==0,1))) = [];
        
        SINR(k,s) = p*eta(k)*bb'*(G\bb);
    end
end

SE = log2(1+SINR(:,:));











