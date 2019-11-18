function [Me] = Cal_metrics(Re,Unmix_s,sOri,SetupStruc)
sign_plot = 1;%choose if plot with results
sign_onLine = 0;
%%%%%%%%%%%%%%%
Num = size(sOri,2);
InSIR = zeros(Num,1);
OutSIR = zeros(Num,1);
for i = 1:Num
    InSIR(i) = CSNR(sOri,i);
end
Me.InSIR = InSIR;
method_num = length(Re.method);
for i = 1:method_num
    method_name = Re.method{i};
    eval(strcat('ReStruc = Re.',method_name,';'));
    if(isfield(ReStruc,'S'))
        eval(strcat('[SDR,SIR,SAR,perm]=bss_eval_sources(Re.',method_name,'.S,sOri);'));
        Us = FrePro(Unmix_s,eval(strcat('Re.',method_name,'.W;')),eval(strcat('SetupStruc.',method_name,';')));
        for j = 1:Num
            OutSIR(j) = CSNR(Us(:,:,j),j);
        end
        ImproSIR = OutSIR-InSIR;
        eval(strcat('Me.',method_name,'_ImproSIR=ImproSIR;'));
        eval(strcat('Me.',method_name,'_OutSIR=OutSIR;'));
        eval(strcat('Me.',method_name,'_SDR=SDR;'));
        eval(strcat('Me.',method_name,'_SIR=SIR;'));
        eval(strcat('Me.',method_name,'_SAR=SAR;'));
        if(sign_plot==1)
            S_c = zeros(size(Us,1),Num);
            for j = 1:Num
                S_c(:,j) = sum(Us(:,:,j),2);
            end
            autoPlot(S_c,method_name,SetupStruc.fs,roundn([ImproSIR OutSIR SDR SIR SAR],-2));
        end
    end
    if(isfield(ReStruc,'OnS'))
        S = ReStruc.OnS;
        [oSDR,oSIR,oSAR,~]=bss_eval_sources(S,sOri);
        eval(strcat('Me.',method_name,'_oSDR=oSDR;'));
        eval(strcat('Me.',method_name,'_oSIR=oSIR;'));
        eval(strcat('Me.',method_name,'_oSAR=oSAR;'));
        if(sign_onLine==1)
            eval(strcat('Struc = SetupStruc.',method_name,';'));
            K = Struc.K;
            hop = Struc.hop;
            L = Struc.L;
            Len = (L-1)*hop+K;
            for j = 1:size(sOri,2)
                Sf(j,:,:) = enframe(S(:,j),Len,hop)';
                So(j,:,:) = enframe(sOri(:,j),Len,hop)';
            end
            frame_N = size(Sf,3);
            Num = size(sOri,2);
            SDRf = zeros(Num,frame_N);
            SIRf = zeros(Num,frame_N);
            SARf = zeros(Num,frame_N);
            for j = 1:frame_N
                Sf_temp = Sf(:,:,j)';
                So_temp = So(:,:,j)';
                [SDRf(:,j),SIRf(:,j),SARf(:,j),~]=bss_eval_sources(Sf_temp,So_temp);
            end
            autoPlot([mean(SDRf);mean(SIRf);mean(SARf)]',strcat(method_name,'_online'),SetupStruc.fs/(Len-hop),...
                roundn(mean([mean(SDRf);mean(SIRf);mean(SARf)],2),-2));
        end
        if(sign_plot==1)
            autoPlot(S,strcat(method_name,'_onLine'),SetupStruc.fs,roundn([oSDR oSIR oSAR],-2));
        end
    end   
end

return