function [Re,SetupStruc] = OnProcess(method_in,s,SetupStruc)
%%%%%%% This function is a transit of method selecting which calculates transfer function, and supply the interface of the method selection
%%%% Process 'method_in' to a row method(cell)
method = {};
cell_num = size(method_in,1);
if(ismember({'ALL'},method_in{cell_num}) && method_in{cell_num*2}==1)
    method ={'ALL'};
else 
    for i = 1:cell_num
        if(ismember({'ALL'},method_in{i}) && method_in{cell_num+i}==1)
            method ={'ALL'};
            break;
        else
            if(method_in{cell_num+i}==1)
                method = [method;method_in{i}];
            end
        end
    end
end
%%%%
Re.method = {};
sign_plotOnline = 0;      % Choosing whether plot the online processing results
sign_plot = 0;      % Choosing whether plot the offline processing results
sign_write = 1;      % Choosing whether generate the .wav file of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('DSB',method) || ismember('ALL',method))
    SetupStruc.DSB.K = 512;                            %%%%% Adjustable coeficient 'K'
    SetupStruc.DSB.hop = round(SetupStruc.DSB.K/2);        %%%%% Adjustable coeficient 'hop'
    SetupStruc.DSB.L = 1;
    Transfer = Cal_transfer(SetupStruc,'DSB');
    [Re.DSB.S,Re.DSB.W,SetupStruc] = Process_DSB(s,Transfer,SetupStruc);
    [Re.DSB.OnS,SetupStruc.DSB] = SimuOnline(s,Transfer,SetupStruc.DSB,'DSB');
    Re.method = [Re.method;'DSB'];
    if(sign_plotOnline == 1)
        autoPlot(Re.DSB.OnS,'DSB_online',SetupStruc.fs);
    end
    if(sign_plot == 1)
        autoPlot(Re.DSB.S,'DSB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.DSB.S,SetupStruc.fs,'DSB');
        autoWrite(Re.DSB.OnS,SetupStruc.fs,'DSB_online');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('DSB_Mask',method) || ismember('ALL',method))
    SetupStruc.DSB_Mask.K = 512;                            
    SetupStruc.DSB_Mask.hop = round(SetupStruc.DSB_Mask.K/2);        
    SetupStruc.DSB_Mask.L = 1;
    Transfer = Cal_transfer(SetupStruc,'DSB_Mask');
%     [Re.DSB_Mask.S,Re.DSB_Mask.W,SetupStruc] = Process_DSB_Mask(s,Transfer,SetupStruc);
    [Re.DSB_Mask.OnS,SetupStruc.DSB_Mask] = SimuOnline(s,Transfer,SetupStruc.DSB_Mask,'DSB_Mask');
    Re.method = [Re.method;'DSB_Mask'];
    if(sign_plotOnline == 1)
        autoPlot(Re.DSB_Mask.OnS,'DSB_Mask',SetupStruc.fs);
    end
    if(sign_plot == 1)
        autoPlot(Re.DSB_Mask.S,'DSB_Mask',SetupStruc.fs);
    end
    if(sign_write == 1)
%         autoWrite(Re.DSB_Mask.S,SetupStruc.fs,'DSB_Mask');
        autoWrite(Re.DSB_Mask.OnS,SetupStruc.fs,'DSB_Mask');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR',method) || ismember('ALL',method))
    SetupStruc.MVDR.K = 512;                            
    SetupStruc.MVDR.hop = round(SetupStruc.MVDR.K/2);        
%     T = 1;
%     SetupStruc.MVDR.InL = ceil(T*SetupStruc.fs/SetupStruc.MVDR.hop);
    SetupStruc.MVDR.L = 1;
%     SetupStruc.MVDR.alpha = SetupStruc.MVDR.InL/(SetupStruc.MVDR.InL+SetupStruc.MVDR.L);
    SetupStruc.MVDR.alpha = 1;
    Transfer = Cal_transfer(SetupStruc,'MVDR');
%     [Re.MVDR.S,Re.MVDR.W,SetupStruc] = Process_MVDR_Search(s,Transfer,SetupStruc);
    [Re.MVDR.OnS,SetupStruc.MVDR] = SimuOnline(s,Transfer,SetupStruc.MVDR,'MVDR');
    Re.method = [Re.method;'MVDR'];
    if(sign_plotOnline == 1)
        autoPlot(Re.MVDR.OnS,'MVDR_online',SetupStruc.fs);
    end
    if(sign_plot == 1)
        autoPlot(Re.MVDR.S,'MVDR',SetupStruc.fs);
    end
    if(sign_write == 1)
%         autoWrite(Re.MVDR.S,SetupStruc.fs,'MVDR');
        autoWrite(Re.MVDR.OnS,SetupStruc.fs,'MVDR_online');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('maxSNR',method) || ismember('ALL',method))
    SetupStruc.maxSNR.K = 2048;                            
    SetupStruc.maxSNR.hop = round(SetupStruc.maxSNR.K/4);        
    SetupStruc.maxSNR.L = 1;
    Transfer = Cal_transfer(SetupStruc,'maxSNR');
    [Re.maxSNR.S,Re.maxSNR.W,SetupStruc] = Process_maxSNR(s,Transfer,SetupStruc);
    if(isfield(SetupStruc,'unS'))
        SetupStruc.maxSNR.unS = SetupStruc.unS;
    end
    [Re.maxSNR.OnS,SetupStruc.maxSNR] = SimuOnline(s,Transfer,SetupStruc.maxSNR,'maxSNR');
    Re.method = [Re.method;'maxSNR'];
    if(sign_plotOnline == 1)
        autoPlot(Re.maxSNR.OnS,'maxSNR',SetupStruc.fs);
    end
    if(sign_plot == 1)
        autoPlot(Re.maxSNR.S,'maxSNR',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.maxSNR.S,SetupStruc.fs,'maxSNR');
        autoWrite(Re.maxSNR.OnS,SetupStruc.fs,'maxSNR_online');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('LCMV',method) || ismember('ALL',method))
    SetupStruc.LCMV.K = 512;                            
    SetupStruc.LCMV.hop = round(SetupStruc.LCMV.K/2);        
%     T = 0.1;
%     SetupStruc.LCMV.InL = ceil(T*SetupStruc.fs/SetupStruc.LCMV.hop);
    SetupStruc.LCMV.L = 1;
%     SetupStruc.LCMV.alpha = SetupStruc.LCMV.InL/(SetupStruc.LCMV.InL+SetupStruc.LCMV.L);
    SetupStruc.LCMV.alpha = 1;
    Transfer = Cal_transfer(SetupStruc,'LCMV');
    [Re.LCMV.S,Re.LCMV.W,SetupStruc] = Process_LCMV(s,Transfer,SetupStruc);
    [Re.LCMV.OnS,SetupStruc.LCMV] = SimuOnline(s,Transfer,SetupStruc.LCMV,'LCMV');
    Re.method = [Re.method;'LCMV'];
    if(sign_plotOnline == 1)
        autoPlot(Re.LCMV.OnS,'LCMV_online',SetupStruc.fs);
    end
    if(sign_plot == 1)
        autoPlot(Re.LCMV.S,'LCMV',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.LCMV.S,SetupStruc.fs,'LCMV');
        autoWrite(Re.LCMV.OnS,SetupStruc.fs,'LCMV_online');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('ICA',method) || ismember('ALL',method))
    SetupStruc.ICA.K = 1024;                            
    SetupStruc.ICA.hop = round(SetupStruc.ICA.K/4);        
    Transfer = Cal_transfer(SetupStruc,'ICA');
    [Re.ICA.S,Re.ICA.W,SetupStruc] = Process_ICA(s,Transfer,SetupStruc);
    Re.method = [Re.method;'ICA'];
    if(sign_plot == 1)
        autoPlot(Re.ICA.S,'ICA',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.ICA.S,SetupStruc.fs,'ICA');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;