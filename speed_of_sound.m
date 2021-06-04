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
%% ファイル名
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
%% 周波数
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]


%% 実験No
No=str2double(path_name1(1:2));

%% マイクの位置
pos=str2double(path_name2(8:9));
pos=rem(pos,16);
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
v=zeros(1,16);

%% 温度
if No<9
    No=1;
elseif 8 < No&&No <17
    No=9;
elseif 16 < No&&No <25
    No=17;
else
    No=25;
end
T=tem((S2(No).x____testo____________+S2(No).Var11)/2);

%% グラフデータ
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
    
    
    
    % 時間窓設定
    idxmbase=fix((leng(S1(pos).xpos,S1(pos).ypos,S1(ix).xpos,S1(ix).ypos)*192000/T-50+0.18*192))+1;
    if idxmbase < 0
        idxmbase = 1;
    end
    time.spos(1,ix)=idxmbase;
    time.s(1,ix)=xtime(idxmbase);
    time.e(1,ix)=xtime(idxmbase+350);
    
    % 振幅最大値
    [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+349)));
    tpos=idxmbase+tpos;
    tmax.tpos(1,ix)=tpos;
    tmax.ty(1,ix)=tymax;
    tmax.t(1,ix)=xtime(tpos);
    
    % 伝搬時間
    pt.pos(1,ix)=find(abs(wdm(idxmbase:idxmbase+349))>(tymax*0.05), 1 )+idxmbase;
    pt.t(1,ix)=xtime(pt.pos(1,ix))-0.18;
    
    %音速
    v(1,ix)=speed(leng(S1(pos).xpos,S1(pos).ypos,S1(ix).xpos,S1(ix).ypos),pt.t(1,ix));
    
    
    figure(1);subplot(4,4,ix)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([time.s(1,ix),time.e(1,ix)]);ylim([-ymax ymax]);
    hold on
    xline(xtime(tpos),'-r');
    xline(xtime(pt.pos(1,ix)),'-b');
    hold off
    ylabel('Amplitude')
    
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(time.s(1,ix),1.4,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)

    figure(2);subplot(4,4,ix)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([time.s(1,ix),time.e(1,ix)]);ylim([0 1.2]);
    hold on
    xline(xtime(tpos),'-r');
    xline(xtime(pt.pos(1,ix)),'-b');
    hold off
    ylabel('Amplitude')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(time.s(1,ix),1.1,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
    figure(3);subplot(4,4,ix)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([time.s(1,ix),time.e(1,ix)]);ylim([-3.99 3.99]);
    ylabel('Phase[rad]')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(time.s(1,ix),3.5,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
end


%% templeture
function T=tem(t)
    T=331.5+0.6*t;
end

%% length
function len=leng(x1,y1,x2,y2)
    l=(x1-x2)^2+(y1-y2)^2;
    len=sqrt(l);
end

%% 音速
function v=speed(l,t)
    v=l/t;
end


