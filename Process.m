function [Re,SetupStruc] = Process(method_in,s,SetupStruc)
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
sign_plot = 0;      % Choosing whether plot the offline processing results
sign_write = 0;      % Choosing whether generate the .wav file of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('DSB',method) || ismember('ALL',method))
    SetupStruc.DSB.K = 512;                            %%%%% Adjustable coeficient 'K'
    SetupStruc.DSB.hop = round(SetupStruc.DSB.K/2);        %%%%% Adjustable coeficient 'hop'
    Transfer = Cal_transfer(SetupStruc,'DSB');
    [Re.DSB.S,Re.DSB.W,SetupStruc] = Process_DSB(s,Transfer,SetupStruc);
    Re.method = [Re.method;'DSB'];
    if(sign_plot == 1)
        autoPlot(Re.DSB.S,'DSB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.DSB.S,SetupStruc.fs,'DSB');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('DSB_Mask',method) || ismember('ALL',method))
    SetupStruc.DSB_Mask.K = 512;                           
    SetupStruc.DSB_Mask.hop = round(SetupStruc.DSB_Mask.K/2);        
    Transfer = Cal_transfer(SetupStruc,'DSB_Mask');
    [Re.DSB_Mask.S,Re.DSB_Mask.W,SetupStruc] = Process_DSB_Mask(s,Transfer,SetupStruc);
    Re.method = [Re.method;'DSB_Mask'];
    if(sign_plot == 1)
        autoPlot(Re.DSB_Mask.S,'DSB_Mask',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.DSB_Mask.S,SetupStruc.fs,'DSB_Mask');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR',method) || ismember('ALL',method))
    SetupStruc.MVDR.K = 512;                           
    SetupStruc.MVDR.hop = round(SetupStruc.MVDR.K/2);        
    Transfer = Cal_transfer(SetupStruc,'MVDR');
    [Re.MVDR.S,Re.MVDR.W,SetupStruc] = Process_MVDR(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR.S,'MVDR',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR.S,SetupStruc.fs,'MVDR');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR_ESB',method) || ismember('ALL',method))
    SetupStruc.MVDR_ESB.K = 512;                            
    SetupStruc.MVDR_ESB.hop = round(SetupStruc.MVDR_ESB.K/2);        
    Transfer = Cal_transfer(SetupStruc,'MVDR_ESB');
    [Re.MVDR_ESB.S,Re.MVDR_ESB.W,SetupStruc] = Process_MVDR_ESB(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR_ESB'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR_ESB.S,'MVDR_ESB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR_ESB.S,SetupStruc.fs,'MVDR_ESB');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR_AESB',method) || ismember('ALL',method))
    SetupStruc.MVDR_AESB.K = 512;                            
    SetupStruc.MVDR_AESB.hop = round(SetupStruc.MVDR_AESB.K/2);       
    Transfer = Cal_transfer(SetupStruc,'MVDR_AESB');
    [Re.MVDR_AESB.S,Re.MVDR_AESB.W,SetupStruc] = Process_MVDR_AESB(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR_AESB'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR_AESB.S,'MVDR_AESB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR_AESB.S,SetupStruc.fs,'MVDR_AESB');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR_PCA',method) || ismember('ALL',method))
    SetupStruc.MVDR_PCA.K = 512;                            
    SetupStruc.MVDR_PCA.hop = round(SetupStruc.MVDR_PCA.K/2);        
    Transfer = Cal_transfer(SetupStruc,'MVDR_PCA');
    [Re.MVDR_PCA.S,Re.MVDR_PCA.W,SetupStruc] = Process_MVDR_PCA(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR_PCA'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR_PCA.S,'MVDR_PCA',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR_PCA.S,SetupStruc.fs,'MVDR_PCA');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR_Search',method) || ismember('ALL',method))
    SetupStruc.MVDR_Search.K = 512;                            
    SetupStruc.MVDR_Search.hop = round(SetupStruc.MVDR_Search.K/2);        
    Transfer = Cal_transfer(SetupStruc,'MVDR_Search');
    [Re.MVDR_Search.S,Re.MVDR_Search.W,SetupStruc] = Process_MVDR_Search(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR_Search'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR_Search.S,'MVDR_Search',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR_Search.S,SetupStruc.fs,'MVDR_Search');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('MVDR_AESB_Search',method) || ismember('ALL',method))
    SetupStruc.MVDR_AESB_Search.K = 512;                            
    SetupStruc.MVDR_AESB_Search.hop = round(SetupStruc.MVDR_AESB_Search.K/2);        
    Transfer = Cal_transfer(SetupStruc,'MVDR_AESB_Search');
    [Re.MVDR_AESB_Search.S,Re.MVDR_AESB_Search.W,SetupStruc] = Process_MVDR_AESB_Search(s,Transfer,SetupStruc);
    Re.method = [Re.method;'MVDR_AESB_Search'];
    if(sign_plot == 1)
        autoPlot(Re.MVDR_AESB_Search.S,'MVDR_AESB_Search',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.MVDR_AESB_Search.S,SetupStruc.fs,'MVDR_AESB_Search');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('LCMV',method) || ismember('ALL',method))
    SetupStruc.LCMV.K = 512;                            
    SetupStruc.LCMV.hop = round(SetupStruc.LCMV.K/2);        
    Transfer = Cal_transfer(SetupStruc,'LCMV');
    [Re.LCMV.S,Re.LCMV.W,SetupStruc] = Process_LCMV(s,Transfer,SetupStruc);
    Re.method = [Re.method;'LCMV'];
    if(sign_plot == 1)
        autoPlot(Re.LCMV.S,'LCMV',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.LCMV.S,SetupStruc.fs,'LCMV');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('LCMV_ESB',method) || ismember('ALL',method))
    SetupStruc.LCMV_ESB.K = 512;                            
    SetupStruc.LCMV_ESB.hop = round(SetupStruc.LCMV_ESB.K/2);        
    Transfer = Cal_transfer(SetupStruc,'LCMV_ESB');
    [Re.LCMV_ESB.S,Re.LCMV_ESB.W,SetupStruc] = Process_LCMV_ESB(s,Transfer,SetupStruc);
    Re.method = [Re.method;'LCMV_ESB'];
    if(sign_plot == 1)
        autoPlot(Re.LCMV_ESB.S,'LCMV_ESB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.LCMV_ESB.S,SetupStruc.fs,'LCMV_ESB');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('LCMV_AESB',method) || ismember('ALL',method))
    SetupStruc.LCMV_AESB.K = 512;                            
    SetupStruc.LCMV_AESB.hop = round(SetupStruc.LCMV_AESB.K/2);        
    Transfer = Cal_transfer(SetupStruc,'LCMV_AESB');
    [Re.LCMV_AESB.S,Re.LCMV_AESB.W,SetupStruc] = Process_LCMV_AESB(s,Transfer,SetupStruc);
    Re.method = [Re.method;'LCMV_AESB'];
    if(sign_plot == 1)
        autoPlot(Re.LCMV_AESB.S,'LCMV_AESB',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.LCMV_AESB.S,SetupStruc.fs,'LCMV_AESB');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('LCMV_Search',method) || ismember('ALL',method))
    SetupStruc.LCMV_Search.K = 512;                            
    SetupStruc.LCMV_Search.hop = round(SetupStruc.LCMV_Search.K/2);        
    Transfer = Cal_transfer(SetupStruc,'LCMV_Search');
    [Re.LCMV_Search.S,Re.LCMV_Search.W,SetupStruc] = Process_LCMV_Search(s,Transfer,SetupStruc);
    Re.method = [Re.method;'LCMV_Search'];
    if(sign_plot == 1)
        autoPlot(Re.LCMV_Search.S,'LCMV_Search',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.LCMV_Search.S,SetupStruc.fs,'LCMV_Search');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('ICA_funda',method) || ismember('ALL',method))
    SetupStruc.ICA_funda.K = 2048;                            
    SetupStruc.ICA_funda.hop = round(SetupStruc.ICA_funda.K/4);        
    Transfer = Cal_transfer(SetupStruc,'ICA_funda');
    [Re.ICA_funda.S,Re.ICA_funda.W,SetupStruc] = Process_ICA_funda(s,Transfer,SetupStruc);
    Re.method = [Re.method;'ICA_funda'];
    if(sign_plot == 1)
        autoPlot(Re.ICA_funda.S,'ICA_funda',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.ICA_funda.S,SetupStruc.fs,'ICA_funda');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('ICA_initial',method) || ismember('ALL',method))
    SetupStruc.ICA_initial.K = 2048;                            
    SetupStruc.ICA_initial.hop = round(SetupStruc.ICA_initial.K/4);        
    Transfer = Cal_transfer(SetupStruc,'ICA_initial');
    [Re.ICA_initial.S,Re.ICA_initial.W,SetupStruc] = Process_ICA_initial(s,Transfer,SetupStruc);
    Re.method = [Re.method;'ICA_initial'];
    if(sign_plot == 1)
        autoPlot(Re.ICA_initial.S,'ICA_initial',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.ICA_initial.S,SetupStruc.fs,'ICA_initial');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('ICA_Sawada',method) || ismember('ALL',method))
    SetupStruc.ICA_Sawada.K = 2048;                            
    SetupStruc.ICA_Sawada.hop = round(SetupStruc.ICA_Sawada.K/4);        
    [Transfer,SetupStruc] = Cal_transfer(SetupStruc,'ICA_Sawada');
    [Re.ICA_Sawada.S,Re.ICA_Sawada.W,SetupStruc] = Process_ICA_Sawada(s,Transfer,SetupStruc);
    Re.method = [Re.method;'ICA_Sawada'];
    if(sign_plot == 1)
        autoPlot(Re.ICA_Sawada.S,'ICA_Sawada',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.ICA_Sawada.S,SetupStruc.fs,'ICA_Sawada');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('IVA',method) || ismember('ALL',method))
    SetupStruc.IVA.K = 2048;                            
    SetupStruc.IVA.hop = round(SetupStruc.IVA.K/4);        
    Transfer = Cal_transfer(SetupStruc,'IVA');
    [Re.IVA.S,Re.IVA.W,SetupStruc] = Process_IVA(s,Transfer,SetupStruc);
    Re.method = [Re.method;'IVA'];
    if(sign_plot == 1)
        autoPlot(Re.IVA.S,'IVA',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.IVA.S,SetupStruc.fs,'IVA');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('AuxIVA',method) || ismember('ALL',method))
    SetupStruc.AuxIVA.K = 2048;                            
    SetupStruc.AuxIVA.hop = round(SetupStruc.AuxIVA.K/4);        
    Transfer = Cal_transfer(SetupStruc,'AuxIVA');
    [Re.AuxIVA.S,Re.AuxIVA.W,SetupStruc] = Process_AuxIVA(s,Transfer,SetupStruc);
    Re.method = [Re.method;'AuxIVA'];
    if(sign_plot == 1)
        autoPlot(Re.AuxIVA.S,'AuxIVA',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.AuxIVA.S,SetupStruc.fs,'AuxIVA');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('OverIVA',method) || ismember('ALL',method))
    SetupStruc.OverIVA.K = 2048;                            
    SetupStruc.OverIVA.hop = round(SetupStruc.OverIVA.K/4);        
    Transfer = Cal_transfer(SetupStruc,'OverIVA');
    [Re.OverIVA.S,Re.OverIVA.W,SetupStruc] = Process_OverIVA(s,Transfer,SetupStruc);
    Re.method = [Re.method;'OverIVA'];
    if(sign_plot == 1)
        autoPlot(Re.OverIVA.S,'OverIVA',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.OverIVA.S,SetupStruc.fs,'OverIVA');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('maxSNR',method) || ismember('ALL',method))
    SetupStruc.maxSNR.K = 4096;                            
    SetupStruc.maxSNR.hop = round(SetupStruc.maxSNR.K/4);        
    Transfer = Cal_transfer(SetupStruc,'maxSNR');
    [Re.maxSNR.S,Re.maxSNR.W,SetupStruc] = Process_maxSNR(s,Transfer,SetupStruc);
    Re.method = [Re.method;'maxSNR'];
    if(sign_plot == 1)
        autoPlot(Re.maxSNR.S,'maxSNR',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.maxSNR.S,SetupStruc.fs,'maxSNR');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dereverberation method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('WPE',method) || ismember('ALL',method))
    SetupStruc.WPE.K = 1024;                            %%%%% Adjustable coeficient 'K'
    SetupStruc.WPE.hop = round(SetupStruc.WPE.K/2);        %%%%% Adjustable coeficient 'hop'
%     Transfer = Cal_transfer(SetupStruc,'WPE');
    [Re.WPE.S,SetupStruc] = Process_WPE(s,SetupStruc);
    Re.method = [Re.method;'WPE'];
    if(sign_plot == 1)
        autoPlot(Re.WPE.S(:,1),'WPE',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.WPE.S,SetupStruc.fs,'WPE');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compound method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ismember('cGMM_maxSNR',method) || ismember('ALL',method))
    SetupStruc.cGMM_maxSNR.K = 512;                            
    SetupStruc.cGMM_maxSNR.hop = round(SetupStruc.cGMM_maxSNR.K/4);        
    Transfer = Cal_transfer(SetupStruc,'cGMM_maxSNR');
    [Re.cGMM_maxSNR.S,Re.cGMM_maxSNR.W,SetupStruc] = Process_cGMM_maxSNR(s,Transfer,SetupStruc);
    Re.method = [Re.method;'cGMM_maxSNR'];
    if(sign_plot == 1)
        autoPlot(Re.cGMM_maxSNR.S,'cGMM_maxSNR',SetupStruc.fs);
    end
    if(sign_write == 1)
        autoWrite(Re.cGMM_maxSNR.S,SetupStruc.fs,'cGMM_maxSNR');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return;