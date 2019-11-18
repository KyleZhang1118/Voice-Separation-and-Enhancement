function covMa = cal_covMa(S,K,hop,win)
num = size(S,3);
N = size(S,2);
K_m = K/2+1;
for i = 1:num
    for j = 1:N
        X(:,:,j,i) = fft(enframe(S(:,j,i),win,hop)');
    end
end
covMa = zeros(N,N,K_m,num);
for i = 1:num
    for j = 1:K_m
        temp = permute(X(j,:,:,i),[3 2 1 4]);
        covMa(:,:,j,i) = temp*temp'/size(X,2);
    end
end