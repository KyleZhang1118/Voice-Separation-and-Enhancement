function [Y,W,SetupStruc] = Process_LCMV_Search(s,Transfer,SetupStruc)
K = SetupStruc.LCMV_Search.K;
hop = SetupStruc.LCMV_Search.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.LCMV_Search.win = win;  % Preserve 'win' in 'SetupStruc'
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
theta = 10^-7;
W = zeros(Num,N,K_m);
Num_MVDR = 0;
Num_LCMV = 0;
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LCMV processing, imposing the linear constraints to the angle of sources
    W_f = zeros(N,Num);
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    for j = 1:Num       
        sign_MVDR = 0;
        R = X_f*X_f'/frame_N;
        if rcond(R)<theta
            R = R+eye(N)*min([min(diag(R)) theta]);
        end
        if(rcond(Steer'*Steer)<theta)
            h = Steer(:,j);
            w = R\h/(h'/R*h);
            sign_MVDR = 1;Num_MVDR =Num_MVDR+1;
            WNG_set = 3;
        else
            w = R\Steer/(Steer'/R*Steer);
            w = w(:,j);
            Steer_incov = inv(Steer'*Steer);
            Num_LCMV = Num_LCMV+1;
            WNG_set =8.451;%%% 7 channels WNG limit 10*log10(7), dont using 'R' here 
            WNG_set = min([WNG_set -10*log10(Steer_incov(j,j))]);
        end
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
                R = R-e*eye(N);
                e = e/10;
                signM = 0;
            else
                if(signM==1)
                    e = e*2;
                end
                R = R+e*eye(N);
            end
            if(sign_MVDR==1)
                w = R\h/(h'/R*h);
            else
                R = eye(size(R,1));sign=0;   %%%%%%%%%% Fixed beamforming
                w = R\Steer/(Steer'/R*Steer);
                w = w(:,j);
            end
            wng = -10*log10(diag(w'*w));
        end
        W_f(:,j) = w;
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
% autoPlot(WNG,'LCMV_Search');
return;