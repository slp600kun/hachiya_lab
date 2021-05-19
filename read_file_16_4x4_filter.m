%
% wave ファイルのプロット
%
%%
clearvars
close all

% 図の基本設定
ax.LineWidth=1.2;ax.FontSize=12;
ps.Color='black';ps.LineWidth=1.5;
tx.Interpreter='none';

%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
Fc=str2num(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]

ch_str = num2str(1,'%02d') ;
file(3:4) = ch_str ;
[wdata,Fs] = audioread([path file]) ;

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

xmin=0; xmax=15;% 0-15ms
ymin=0; ymax=1.1;

f1=figure(1);
set(f1,'Position', [700 500 1200 800])

f2=figure(2);
set(f2,'Position', [600 300 1200 800])

f3=figure(3);
set(f3,'Position', [500 100 1200 800])

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
    
    figure(1);subplot(4,4,ix)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([xmin xmax]);ylim([-ymax ymax]);
    ylabel('Amplitude')
    
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
        tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
        xlabel('Time(ms)')
    end
    
    text(10,0.85,['CH : ' ch_str],'FontSize',12)
    
    set(gca,ax)
    
    figure(2);subplot(4,4,ix)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([xmin xmax]);ylim([0 ymax]);
    ylabel('Amplitude')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
        tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
        xlabel('Time(ms)')
    end
    text(10,0.95,['CH : ' ch_str],'FontSize',12)
    set(gca,ax)
    
    figure(3);subplot(4,4,ix)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([xmin xmax]);ylim([-3.99 3.99]);
    ylabel('Phase[rad]')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
        tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
        xlabel('Time(ms)')
    end
    text(10,3.5,['CH : ' ch_str],'FontSize',12)
    set(gca,ax)
    
    
end