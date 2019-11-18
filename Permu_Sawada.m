function P = Permu_Sawada(W,Y,SetupStruc,method)
%%
thA = 20; %% Angle/degree
thU = 3; %% SIR/dB
thU_mul = 10^(thU/10);
interval = 3; %% Field of cor/+-bins
thcor = 0.7; %% Correlation
thcor_single = 0.8;
single_sign = 0; %%select if could judge permutation by one cal of cor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Num = size(W,1);
N = size(W,2);
K_m = size(W,3);
P = zeros(Num,Num,K_m);
Permu = zeros(Num,K_m-1);
Permu_sign = zeros(K_m-1,1);   %%% 1 repersents bin aligned,0 not aligned, -1 cant aliged,2 remains to fresh domian,3 caled with hamonic
Permu_cir = zeros(factorial(Num),K_m-1);
Permu_Num = zeros(K_m-1,1);
Permu_NumL = ones(K_m-1,1)*interval*2;
for i = 1:interval
    Permu_NumL(i) = i+interval-1;
    Permu_NumL(K_m-i) = i+interval-1;
end
Rate_meth = zeros(K_m,4);
%%%% I used the information of angle which should calculated by the inverse of W, this process is time-consuming. 
Angle = SetupStruc.Angle;
Transfer = SetupStruc.Transfer;
if(~exist('Angle','var'))
    sprintf('Here should calculate Angle by W+.');
end
if(~exist('method','var'))
    method = 'all';
end
if(strcmp(method,'DOA'))    
    thA = 360;
    thU_mul = 0;
elseif(strcmp(method,'cor'))
    thcor = -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOA
%%%%This is a simple realization without angle interval detection. 
for i = 2:K_m
    trans = permute(Transfer(i,:,:),[2 3 1]);
    Pat_count = zeros(Num,1);
    for j = 1:Num
        Beam_pat = W(j,:,i)*trans;
        Beam_pat = abs(Beam_pat);
        [max1,pos1] = max(Beam_pat);
        Beam_pat(pos1) = -1;
        max2 = max(Beam_pat);
        if(max2*thU_mul<=max1)
            Permu(j,i-1) = pos1;
            Pat_count(pos1) = Pat_count(pos1)+1;
        end
    end
    if(~isempty(find(Pat_count>1, 1)))  %%%%Judge if there're more than 1 row of W target the same angle.
        for j = 1:Num
            if(Permu(j,i-1)>0)
                if(Pat_count(Permu(j,i-1))>1)
                    Pat_count(Permu(j,i-1)) = 0;
                    Permu(j,i-1) = 0;
                end
            end
        end
    end
    Rate_meth(i-1,1) = length(find(Pat_count==1))/Num;
    zer = find(Pat_count==0);
    if(isempty(zer))
        Permu_sign(i-1) = 2;
    elseif(length(zer)==1)
        Permu(Permu(:,i-1)==0,i-1) = zer;
        Permu_sign(i-1) = 2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation method
if(~strcmp(method,'DOA'))
    if(isempty(find(Permu_sign==2,1)))
        Permu(:,round(K_m/4)) = 1:Num;
        Permu_sign(round(K_m/4)) = 2;
    end
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
end
Permu_pending = find(Permu_sign==2);
while(~isempty(Permu_pending))
    for i = 1:length(Permu_pending)
        x_l = Permu_pending(i)-interval;
        x_h = Permu_pending(i)+interval;
        if(x_l<1)
            x_l = 1;
        end
        if(x_h>K_m-1)
            x_h = K_m-1;
        end
        for j = x_l:x_h
            if(Permu_sign(j)==0)
                cor_cir = zeros(size(Table_Permu,1),1);
                for k = 1:size(Table_Permu,1)
                    if(~isempty(find(Permu(:,j)~=0,1)))
                        verif = find(Permu(:,j)~=0);
                        sign_verif = 0;
                        for k_i = 1:length(verif)
                            if(Table_Permu(k,verif(k_i))~=Permu(verif(k_i),j))
                                sign_verif = 1;
                                break;
                            end
                        end
                        if(sign_verif==1)
                            continue;
                        end
                    end
                    for k_i = 1:Num
                        y1_temp = abs(Y(:,k_i,j+1));
                        index = Permu(:,Permu_pending(i))==Table_Permu(k,k_i);
                        y2_temp = abs(Y(:,index,Permu_pending(i)+1));
                        cor = corr(y1_temp,y2_temp);
                        cor_cir(k) = cor_cir(k)+abs(cor);                       
                    end
                end
%                 cor_cir = cal_cor(Y(:,:,j+1),Y(:,:,Permu_pending(i)+1)); %%%%%%%%%%%%%%
                if(single_sign==1 && max(cor_cir)>=thcor_single*Num)
                    [~,index] = max(cor_cir);
                    Permu(:,j) = Table_Permu(index,:)';   
                    Permu_sign(j) = 2;
                    Rate_meth(j,2) = 0.5;
                else
                    Permu_cir(:,j) = Permu_cir(:,j)*Permu_Num(j)+cor_cir;
                    Permu_Num(j) = Permu_Num(j)+1;                  
                    if(Permu_Num(j)==Permu_NumL(j))%%%%%%%%
                        [cor_in,index] = max(Permu_cir(:,j));
                        if(strcmp(method,'cor') || cor_in>=thcor*Num*Permu_Num(j))                           
                            Permu(:,j) = Table_Permu(index,:)';   
                            Permu_sign(j) = 2;
                            Rate_meth(j,2) = 1;
                        end
                    end
                    Permu_cir(:,j) = Permu_cir(:,j)/Permu_Num(j);
                end
            end
        end
        Permu_sign(Permu_pending(i)) = 1;
%         Permu_cir(:,Permu_pending(i)) = zeros(size(Permu_cir,1),1);
    end
    Permu_pending = find(Permu_sign==2);
    if(isempty(Permu_pending))
        if(find(Permu_sign==0,1))
            cir0 = find(Permu_sign==0);
            icor_max = max(Permu_cir(:,cir0));
            [cor_max,i_max] = max(icor_max);
            if(strcmp(method,'cor') || cor_max>=thcor*Num)
                [~,index] = max(Permu_cir(:,cir0(i_max)));
                Permu(:,cir0(i_max)) = Table_Permu(index,:)';
                Permu_sign(cir0(i_max)) = 2;
                Rate_meth(cir0(i_max),2) = 1;
                Permu_pending = cir0(i_max);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Harmonic method
if(strcmp(method,'all'))
    Harm_lim = floor((K_m-2)/3);
    Permu_Harm = zeros(size(Table_Permu,1),Harm_lim);
    Permu_NumH = zeros(Harm_lim,1);
    for i = 2:Harm_lim
        Harmo = [2*i-1 2*i 2*i+1 3*i-1 3*i 3*i+1];
        Harmo_Num = 6;
        if(i==2)
            Harmo = [2*i-1 2*i 2*i+1 3*i 3*i+1];
            Harmo_Num = 5;
        end
        if(Permu_sign(i)==0)
            for j = 1:Harmo_Num
                cor_cir = zeros(size(Table_Permu,1),1);
                if(Permu_sign(Harmo(j))==1)
                    for k = 1:size(Table_Permu,1)
                        if(~isempty(find(Permu(:,i)~=0,1)))
                            verif = find(Permu(:,i)~=0);
                            sign_verif = 0;
                            for k_i = 1:length(verif)
                                if(Table_Permu(k,verif(k_i))~=Permu(verif(k_i),i))
                                    sign_verif = 1;
                                    break;
                                end
                            end
                            if(sign_verif==1)
                                continue;
                            end
                        end
                        for k_i = 1:Num
                            y1_temp = abs(Y(:,k_i,i+1));
                            index = Permu(:,Harmo(j))==Table_Permu(k,k_i);
                            y2_temp = abs(Y(:,index,Harmo(j)+1));
                            cor = corr(y1_temp,y2_temp);
                            cor_cir(k) = cor_cir(k)+abs(cor);
                        end
                    end
                    if(single_sign==1 && max(cor_cir)>=thcor_single*Num)
                        [~,index] = max(cor_cir);
                        Permu(:,i) = Table_Permu(index,:)';   
                        Permu_sign(i) = 2;
                        Rate_meth(i,3) = 0.5;
                    else
                        Permu_Harm(:,i) = Permu_Harm(:,i)+cor_cir;
                        Permu_NumH(i) = Permu_NumH(i)+1;                                         
                    end
                end
            end
            cor_in = max(Permu_cir(:,i));
            [cor_ha,index] = max(Permu_Harm(:,i));
            if(cor_ha/(Permu_NumH(i)+eps)>cor_in)                           
                Permu(:,i) = Table_Permu(index,:)';   
                Permu_sign(i) = 2;
                Rate_meth(i,3) = 1;
            end
        end       
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Recal cor
    Permu_pending = find(Permu_sign==2);
    while(~isempty(Permu_pending))
        for i = 1:length(Permu_pending)
            x_l = Permu_pending(i)-interval;
            x_h = Permu_pending(i)+interval;
            if(x_l<1)
                x_l = 1;
            end
            if(x_h>K_m-1)
                x_h = K_m-1;
            end
            for j = x_l:x_h
                if(Permu_sign(j)==0)
                    cor_cir = zeros(size(Table_Permu,1),1);
                    for k = 1:size(Table_Permu,1)
                        if(~isempty(find(Permu(:,j)~=0,1)))
                            verif = find(Permu(:,j)~=0);
                            sign_verif = 0;
                            for k_i = 1:length(verif)
                                if(Table_Permu(k,verif(k_i))~=Permu(verif(k_i),j))
                                    sign_verif = 1;
                                    break;
                                end
                            end
                            if(sign_verif==1)
                                continue;
                            end
                        end
                        for k_i = 1:Num
                            y1_temp = abs(Y(:,k_i,j+1));
                            index = Permu(:,Permu_pending(i))==Table_Permu(k,k_i);
                            y2_temp = abs(Y(:,index,Permu_pending(i)+1));
                            cor = corr(y1_temp,y2_temp);
                            cor_cir(k) = cor_cir(k)+abs(cor);
                        end
                    end
%                 cor_cir = cal_cor(Y(:,:,j+1),Y(:,:,Permu_pending(i)+1)); %%%%%%%%%%%%%%
                    if(single_sign==1 && max(cor_cir)>=thcor_single*Num)
                        [~,index] = max(cor_cir);
                        Permu(:,j) = Table_Permu(index,:)';   
                        Permu_sign(j) = 2;
                        Rate_meth(j,4) = 0.5;
                    else
                        Permu_cir(:,j) = Permu_cir(:,j)*Permu_Num(j)+cor_cir;
                        Permu_Num(j) = Permu_Num(j)+1;                  
                        if(Permu_Num(j)==Permu_NumL(j))%%%%%%%%
                            [~,index] = max(Permu_cir(:,j));                        
                            Permu(:,j) = Table_Permu(index,:)';   
                            Permu_sign(j) = 2;
                            Rate_meth(j,4) = 1;
                        end
                        Permu_cir(:,j) = Permu_cir(:,j)/Permu_Num(j);
                    end
                end
            end
            Permu_sign(Permu_pending(i)) = 1;
%         Permu_cir(:,Permu_pending(i)) = zeros(size(Permu_cir,1),1);
        end
        Permu_pending = find(Permu_sign==2);
        if(isempty(Permu_pending))
            if(find(Permu_sign==0,1))
                cir0 = find(Permu_sign==0);
                icor_max = max(Permu_cir(:,cir0));
                [~,i_max] = max(icor_max);
                [~,index] = max(Permu_cir(:,cir0(i_max)));
                Permu(:,cir0(i_max)) = Table_Permu(index,:)';
                Permu_sign(cir0(i_max)) = 2;
                Rate_meth(cir0(i_max),4) = 1;
                Permu_pending = cir0(i_max);               
            end
        end
    end
end      
%% Generate P
for i = 2:K_m
    W_inv = pinv(W(:,:,i));
    if(~isempty(find(Permu(:,i-1)==0,1)))
        P_sign = zeros(Num,1);
        for j = 1:Num
            if(Permu(j,i-1)~=0)
                P_sign(Permu(j,i-1)) = 1;
            end
        end
        for j = 1:Num
            if(P_sign(j)==0)
                j_sign = j;
                P_sign(j_sign) = 1;
                break;
            end
        end
        for j = 1:Num
            if(Permu(j,i-1)==0)
                P(j,j_sign,i) = W_inv(1,j_sign);
                for k = j_sign:Num
                    if(P_sign(k)==0)
                        j_sign = k;
                        break;
                    end
                end
            else
                P(j,Permu(j,i-1),i) = W_inv(1,Permu(j,i-1));
            end
        end
    else
        for j = 1:Num
            P(j,Permu(j,i-1),i) = W_inv(1,Permu(j,i-1));
        end
    end
end
%%
autoPlot(Rate_meth,1);                
return;






