function [Y,W,SetupStruc] = Process_MVDR_AESB_Search(s,Transfer,SetupStruc)
K = SetupStruc.MVDR_AESB_Search.K;
hop = SetupStruc.MVDR_AESB_Search.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.MVDR_AESB_Search.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(size(X,1),size(X,2),Num);
WNG =zeros(K_m,Num);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-4;
WNG_set = -10;
W = zeros(Num,N,K_m);
sign_zero=0;sign_N=0;sign_Num=0;
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MVDR processing
    W_f = zeros(N,Num);
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    [E,D] = PCA(X_f);
    Num_PCA = size(D,1);
    [~,adaD] = PCA(X_f,1,N,D(1,1)*theta);
    aNum = size(adaD,1);
    if(aNum>Num) aNum = Num;sign_N=sign_N+1;end
    if(aNum==Num) sign_Num = sign_Num+1; end
    if(Num_PCA<1) 
        R_inv = eye(N);
        for j = 1:Num
            h = Steer(:,j);
            W_f(:,j) = R_inv*h/(h'*R_inv*h);
        end
        sign_zero=sign_zero+1;
    else
        if(aNum>Num_PCA)
            aNum = Num_PCA;
        end
        E = E(:,1:aNum);
        D = D(1:aNum,1:aNum);
%         if(~isempty(find(diag(D)<theta, 1)))
%             D = D+min(diag(D))*eye(aNum);
%         end
%         R_inv = E/D*E';
        for j = 1:Num
            h = Steer(:,j);
            w = E/D*E'*h/(h'*E/D*E'*h);
            wng = -10*log10(diag(w'*w));
            e = 0.05;
            if(WNG_set-wng>0)
                sign = 1;
            else 
                sign = 0;
            end
            signM = 1;
            while(sign==1 && abs(wng-WNG_set)>0.5)
                if(wng>WNG_set)
                    D = D-e*eye(aNum);
                    e = e/10;
                    signM = 0;
                else
                    if(signM==1)
                        e = e*2;
                    end
                    D = D+e*eye(aNum);
                end
                w = E/D*E'*h/(h'*E/D*E'*h);
                wng = -10*log10(diag(w'*w));
            end  
            W_f(:,j) = w;
        end
    end
    W_f = W_f';
    W(:,:,i) = W_f;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WNG(i,:) = -10*log10(diag(W_f*W_f')');
    Y_ = W_f*X_f;
    Y_f(i,:,:) = Y_.';
    if(i~=K_m)
        Y_f(K+2-i,:,:) = Y_';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end
% autoPlot(WNG,'MVDR_AESB_Search');
return;