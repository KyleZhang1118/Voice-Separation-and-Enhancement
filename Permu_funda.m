function P = Permu_funda(W,Y,Y_ori)
Num = size(W,1);
N = size(W,2);
K_m = size(W,3);
P = zeros(Num,Num,K_m);
Permu = zeros(Num,K_m-1);
%%%%% Generate the table of permutation circumstance
Table_Permu = zeros(factorial(Num),Num);
Table_Permu(1,:) = 1:Num;
occup = ones(Num,1);
pin = Num;
for i = 2:factorial(Num)
    while(isempty(find(occup==0,1)) || (pin>0 && Table_Permu(i-1,pin)>find(occup==0, 1, 'last' )))
        occup(Table_Permu(i-1,pin)) = 0;
        pin = pin-1;
    end
    index = find(occup==0);
    Table_Permu(i,pin) = index(find(index>Table_Permu(i-1,pin),1));
    occup(Table_Permu(i-1,pin)) = 0;
    occup(Table_Permu(i,pin)) = 1;
    Table_Permu(i,pin+1:Num) = find(occup==0);
    occup = ones(Num,1);
    if(pin>1)
        Table_Permu(i,1:pin-1) = Table_Permu(i-1,1:pin-1);
    end
    pin = Num;
end
for i = 1:K_m-1
    cor_cir = zeros(size(Table_Permu,1),1);
    for k = 1:size(Table_Permu,1)        
        for k_i = 1:Num
            y1_temp = abs(Y(:,k_i,i+1));
            y2_temp = abs(Y_ori(:,Table_Permu(k,k_i),i+1));
            cor = corr(y1_temp,y2_temp);
            cor_cir(k) = cor_cir(k)+abs(cor);                       
         end
    end
    [~,index] = max(cor_cir);
    Permu(:,i) = Table_Permu(index,:)';
end
for i = 2:K_m
    W_inv = pinv(W(:,:,i));
    for j = 1:Num
        P(j,Permu(j,i-1),i) = W_inv(1,Permu(j,i-1));
    end
end 
return;