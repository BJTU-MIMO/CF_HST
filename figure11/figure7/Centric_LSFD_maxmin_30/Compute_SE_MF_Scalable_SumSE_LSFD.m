function [SE] = Compute_SE_MF_Scalable_SumSE_LSFD(R,Q,hmean,I,ID,K,L,M,p,sigma2,D)

[a_weight] = functionGetWeight(R,Q,hmean,I,ID,K,L,M,p,sigma2,D);

% ---------------------------------------------------
a_kls = a_weight;
%a_kls = 1/L*ones(K,L,M);
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

DS = abs(term1(:,:)).^2; % k,s
IT1 = squeeze(sum(real(term2(:,:,:,:)),4));   % k,i,s
IT2 = squeeze(sum(abs(term3(:,:,:)).^2,3)); % k,s
IT3 = squeeze(sum(abs(term4(:,:,:,:)).^2,4)); %k,i,s
IT4 = sigma2*real(term5(:,:)); % k,s

cc = zeros(K,K,M);
for s = 1:M
    for k = 1:K
        for i = 1:K
            if i == k
                cc(i,k,s) = IT1(k,i,s) + IT2(k,s);
            else
                cc(i,k,s) = IT1(k,i,s) + IT3(k,i,s);
            end
        end
    end
end

% -------------------------------------------------------------------

SE = zeros(K,M);
for s = 1:M
    %Initialize the current objective function
    objec_new = inf;
    %Initialize the power control coefficients randomly
    eta = p*rand(K,1);
    %Prepare arrays to store e_k, u_k, and d_k in Algoirthm 7.2
    eee = zeros(K,M);
    uuu = zeros(K,M);
    ddd = zeros(K,M);
    %Intialize the difference between current and previous objective values
    diff = 100;
    
    %Initizalize the iteration counter to zero
    iterr = 0;
    %Go through the algorithm steps if the objective function is improved
    %more than 0.0001 or not improved at all
    while (diff>0.0001) || (diff<0)
        %disp(['iterr------------------- ' num2str(iterr) ' out of ' num2str(diff)]);
        %Increase iteration counter by one
        iterr = iterr+1;
        %Update the previous objective value by the current objective value
        objec_old = objec_new;
        
        %Go through all UEs
        for k = 1:K
            %Compute the numerator and denominator in Line 4 of Algorithm 7.2
            numm = sqrt(eta(k)*DS(k,s));
            denomm = cc(:,k,s)'*eta(:)+IT4(k,s)+eta(k)*DS(k,s);
            
            %Update u_k according to Line 4
            uuu(k,s) = numm/denomm;
            %Update e_k and d_k as in Line 5
            eee(k,s) = 1-abs(numm)^2/denomm;
            ddd(k,s) = 1/eee(k,s);
            
        end
        
        %Go through all UEs
        for k = 1:K
            %Compute p_k as in Line 6
            numm = DS(k,s)*ddd(k,s)*uuu(k,s)^2;
            denomm = numm;
            for i = 1:K
                denomm = denomm+ddd(i,s)*abs(uuu(i,s))^2*cc(k,i,s);
            end
            
            eta(k) = min(p, ddd(k,s)*numm/(denomm^2));
        end
        
        %Update the current objective value
        objec_new = sum(ddd(:,s).*eee(:,s)-log(ddd(:,s)));
        %Obtain the difference between current and previous objective values
        diff = objec_old - objec_new;
        
    end
    
    SE(:,s) = log2(ddd(:,s));
end

