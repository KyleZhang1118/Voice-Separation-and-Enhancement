function Y = WPE_DM(X)
Kl = 10;
D = 2;
iter_max = 10;
K_m = size(X,1);
frame_N = size(X,2);
channel_N = size(X,3);
frame_N2 = frame_N-Kl-D;
Y_= zeros(Kl*channel_N,frame_N2);
Y_cor = zeros(Kl*channel_N,Kl*channel_N,frame_N2);
% G = zeros(Kl*channel_N,channel_N,K_m-1);
% s_ = zeros(channel_N,frame_N2);
A = zeros(iter_max,K_m);
Y = zeros(size(X));
theta = 10^-4;
Obj = zeros(100,K_m)-1;
for f = 2:K_m
    X_f = permute(X(f,:,:),[3,2,1]);
    g = zeros(Kl*channel_N,channel_N);
    %%%%%%%%%%%%%%%生成Y_长矢量,第一列对应Kl+D+1列X_f
    for i = 1:frame_N2
        for j = 1:Kl
            a = (Kl-j)*channel_N+1;
            Y_(a:a+channel_N-1,i) = X_f(:,j+i-1);
        end
        Y_cor(:,:,i) = Y_(:,i)*Y_(:,i)';
    end
    for i_iter = 1:iter_max
        s_ = X_f(:,Kl+D+1:end)-g'*Y_;
        for i = 1:channel_N
            Y_cor_sum = zeros(Kl*channel_N,Kl*channel_N);
            Y_sum = zeros(Kl*channel_N,1);
            for j = 1:frame_N2
                coe = real(s_(i,j)*s_(i,j)');
                Y_cor_sum = Y_cor_sum+Y_cor(:,:,j)/coe;
                Y_sum = Y_sum+Y_(:,j)*X_f(i,Kl+D+j)'/coe;
                Obj(i_iter,f) = Obj(i_iter,f)+log(coe)+1;
            end
            if(rcond(Y_cor_sum)<theta)
                Y_cor_sum = Y_cor_sum+eye(size(Y_cor_sum))*min(diag(Y_cor_sum));
            end
            g(:,i) = Y_cor_sum\Y_sum;
%             opts.SYM = true;
%             g(:,i) = linsolve(Y_cor_sum,Y_sum,opts);
        end
        Obj(i_iter,f) = Obj(i_iter,f)/frame_N2/channel_N;
        delta_s = s_-(X_f(:,Kl+D+1:end)-g'*Y_);
        A(i_iter,f) = sum(sum(abs(delta_s)))/channel_N/frame_N2;
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s_ = X_f(:,Kl+D+1:end)-g'*Y_;
    Y(f,Kl+D+1:end,:) = s_.';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%补充前几帧的去混响
    Y_head = zeros(Kl*channel_N,Kl+D);
    for i = 1+D:Kl+D
        for j = 1:i-D
            a = (j-1)*channel_N+1;
            Y_head(a:a+channel_N-1,i) = X_f(:,i-D-j+1);
        end
    end
    s_temp = X_f(:,1:Kl+D)-g'*Y_head;
    Y(f,1:Kl+D,:) = s_temp.';
end
return
