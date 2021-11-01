function covMa = cal_covMa4Block(Block,blockWin)
K = size(Block,1);
num = size(Block,4);
N = size(Block,2);
K_m = K/2+1;
X = zeros(size(Block));
for i = 1:num
    X(:,:,:,i) = fft(Block(:,:,i).*blockWin);
end
covMa = zeros(N,N,K_m,num);
for i = 1:num
    for j = 1:K_m
        temp = permute(X(j,:,:,i),[2 3 1 4]);
        covMa(:,:,j,i) = temp*temp'/size(X,3);
    end
end