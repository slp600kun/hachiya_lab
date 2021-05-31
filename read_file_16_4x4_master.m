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
% スプレッドシート
filename1 = '6_sonar_tr4x4_pos_b41.csv';
T1 = readtable(filename1);
S1 = table2struct(T1);
filename2 = '1_meas_item_list_rev.xlsx';
T2 = readtable(filename2);
S2 = table2struct(T2);
%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]
%% マイクの位置
pos=str2double(path(pos_frame+7:pos_frame+8));
pos=rem(str2double(path(pos_frame+7:pos_frame+8)),16);
if pos==0
    pos=16;
end

%%
ch_str = num2str(1,'%02d') ;
file(3:4) = ch_str ;
pos_str = num2str(pos,'%02d');
[wdata,Fs] = audioread([path file]) ;

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

ymin=0; ymax=1.5;
t1=zeros(1,16);
t2=zeros(1,16);
idxm=zeros(1,16);
tpo=zeros(1,16);

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
    
    
    
    % 伝播時間
    idxmbase=fix((leng(S1(pos).xpos,S1(pos).ypos,S1(ix).xpos,S1(ix).ypos)*192/0.338-300+0.18*192))+1;
    if idxmbase < 0
        idxmbase = 1;
    end
    idxm(1,ix)=idxmbase;
    [tymax,tpos]=max(abs(wdm(idxmbase+1:idxmbase+300)));
    tpos=idxmbase+tpos;
    
    tp=xtime(tpos);
    tm=xtime(tpos+300)-0.18;
    t1(1,ix)=tp;
    T1=t1(1,ix);
    t2(1,ix)=tm;
    T2=t2(1,ix);
    
    figure(1);subplot(4,4,ix)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([T1,T2]);ylim([-ymax ymax]);
    hold on
    xline(xtime(tpos),'-r');
    hold off
    ylabel('Amplitude')
    
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,1.4,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)

    figure(2);subplot(4,4,ix)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([T1,T2]);ylim([0 1.2]);
    hold on
    xline(xtime(tpos),'-r');
    hold off
    ylabel('Amplitude')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,1.1,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
    figure(3);subplot(4,4,ix)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([T1,T2]);ylim([-3.99 3.99]);
    ylabel('Phase[rad]')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,3.5,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
end

%% length
function len=leng(x1,y1,x2,y2)
    l=(x1-x2)^2+(y1-y2)^2;
    len=sqrt(l);
end


