function Y = FrePro(X,W,SetupStruc)
K = SetupStruc.K;
hop = SetupStruc.hop;
win = SetupStruc.win;
K_m = K/2+1;
Num = size(X,3);
channel_N = size(X,2);
for i = 1:Num
    for j = 1:channel_N
        S(:,:,j,i) = fft(enframe(X(:,j,i),win,hop)');
    end
end  
frame_N = size(S,2);
Length = (frame_N-1)*hop+K;
Y = zeros(Length,Num,Num);
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num    %%% the cycle process signal by using each 'W'
    for j = 1:Num    %%% the cycle process each source with the same 'W'
        X_temp = zeros(K,frame_N);        
        for f = 2:K_m
            S_f = permute(S(f,:,:,j),[3 2 1 4]);
            X_temp(f,:) = W(i,:,f)*S_f;
            if(f<K_m)
                X_temp(2*K_m-f,:) = conj(X_temp(f,:));
            end
        end
        %%%%%% convert X_temp from frequency-domain to time-domain       
        Y(:,j,i) = overlapadd(real(ifft(X_temp))',win,hop);
    end
end    
return;