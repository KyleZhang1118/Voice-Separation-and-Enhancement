function [Y,W,SetupStruc] = Process_FastICA_HO_Sawada(s,Transfer,SetupStruc)
K = SetupStruc.FastICA_HO_Sawada.K;
hop = SetupStruc.FastICA_HO_Sawada.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.FastICA_HO_Sawada.win = win;  % Preserve 'win' in 'SetupStruc'
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
Y_P = zeros(frame_N,Num,K_m);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-6;
W = zeros(Num,N,K_m);
A = zeros(1001,K/2)-1; %%%% Show the decrease of non-linear correlation, ICA max iterations 1000
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PCA and ICA processing
    [E,D] = PCA(X_f,1,Num);
    V = sqrt(D)\E';
    X_f = V*X_f;
    %%%%%%%% FastICA based on higher-order statistics
    W_f = eye(Num);
    pObj = inf;
    for i_iteration = 1:200
        Obj = 0;
        for i_n = 1:Num
            W_i = W_f(i_n,:);       
            y_ = W_i*X_f;
            y_2 = real(y_).^(2)+imag(y_).^(2);
            %%%%%%%%%%%%%%%%%%%%% There are three different contrast functions.
%             %%%  G1(y) = sqrt(a1+y), g1(y) = 1/(2*sqrt(a1+y))
%             g = 0.5*(0.1+y_2).^(-1/2);
%             g_ = -0.25*(0.1+y_2).^(-3/2);
%             Obj = Obj+sum(sqrt(0.1+y_2))/frame_N;
            %%%  G2(y) = log(a2+y),  g2(y) = 1/(a2+y)
            g = (0.1+y_2).^(-1);
            g_ = -(0.1+y_2).^(-2);
            Obj = Obj+sum(log(0.1+y_2))/frame_N;
%             %%%  G3(y) = 0.5*y^2,    g3(y) = y     (Kurtosis)
%             g = y_2;
%             g_ = ones(1,frame_N);
%             Obj = Obj+0.5*sum(y_2.^2)/frame_N;
            %%%%%%%%%%%%%%%%%%%%%
%             W_i = sum(X_f.*repmat(conj(y_).*g,[Num,1]),2)/frame_N-sum(g+y_2.*g_,2)/frame_N*W_i';
            W_i = (X_f*(conj(y_).*g).'-(sum(g,2)+y_2*g_.')*W_i')/frame_N;
            W_i = W_i/norm(W_i);
            W_f(i_n,:) = W_i';
        end
        W_f = (W_f*W_f')^(-1/2)*W_f;
        dObj = pObj-Obj;
        pObj = Obj;
        A(i_iteration,i-1) = Obj;
        if(abs(dObj)/abs(Obj)<theta)
            break;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Y_,W_ICA,A] = FDICA(X_f,eye(Num),A,i);  %%% 'A', 'i' record the decrease for observation
    Y_ = W_f*X_f;
    W(:,:,i) = W_f*V;
    Y_P(:,:,i) = Y_.';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the ambiguity of permutation and amplitude
P = Permu_Sawada(W,Y_P,SetupStruc,'all');  %%%% Options: 'DOA','cor', 'all'
for i = 2:K_m
    W(:,:,i) = P(:,:,i)*W(:,:,i);
    Y_ = permute(Y_P(:,:,i),[2 1 3]);
    Y_ = P(:,:,i)*Y_;
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
return;