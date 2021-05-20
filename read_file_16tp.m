%
% wave ファイルのプロット
%
%%
clearvars
close all


%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');

pos_freq=strfind(path,'kHz');
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]
pos_frame=strfind(path,'frame');
pos=rem(str2double(path(pos_frame+7:pos_frame+8)),16);
if pos==0
    pos=16;
end

ch_str = num2str(1,'%02d') ;
file(3:4) = ch_str ;
[wdata,Fs] = audioread([path file]) ;

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

xs.tp=zeros(1,16);


%% Filter 設定
iftr=1;

if iftr==1
    Fstop = 7000;
    Fpass = 8000;
    Astop = 30;
    Apass = 0.5;
    
    d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
        'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
        'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');
    
    delay = round(mean(grpdelay(d)))+1;
end

%%
for ix=1:16
    ch_str = num2str(ix,'%02d');
    if ix ~= 1
        file(3:4) = ch_str ;
    end
    
    % Filter 処理
    if iftr==1
        [wdata_tmp,Fs] = audioread([path file]) ;
        wdata=filter(d,wdata_tmp) ;
        wdata(1:end-delay+1) = wdata(delay:end);
    else
        [wdata,Fs] = audioread([path file]);
    end
    
    [b, a] = demod(wdata, Fc, Fs, 'qam');
    wdm = complex(a, b);
    
    
    % xrange
    wdm_abs=abs(wdm);
    for i=1:dnum
        if wdm_abs(i+15,1)-wdm_abs(i,1)>0.04
            x1=xtime(1,i);
            break
        end
    end
    
    % propagation time
    xs.tp(1,ix)=x1;

end

if rem(pos,4)==1
    tp=xs.tp(1,pos+11)/1000;
elseif rem(pos,4)==2
    tp=xs.tp(1,pos+9)/1000;
elseif rem(pos,4)==3
    tp=xs.tp(1,pos+7)/1000;
else
    tp=xs.tp(1,pos+5)/1000;
end