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
T1 = readtable(filename1,'VariableNamingRule','preserve');
S1 = table2struct(T1);
filename2 = '1_meas_item_list_rev.xlsx';
T2 = readtable(filename2,'VariableNamingRule','preserve');
T2.Properties.VariableNames{'・データは2回ずつ取得'} = 'Var1';
T2.Properties.VariableNames{'※湿度はtesto温湿度計を使用した参考値'} = 'Var10';
S2 = table2struct(T2);


%% ファイル名
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
pos_wn=strfind(path,'wn');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);

%% 周波数
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]

%% 実験No
No=str2double(path_name1(1:2));

%% マイクの位置
pos = str2double(path_name2(8:9));
pos = rem(pos,16);
if pos == 0
    pos = 16;
end
oppos=str2double(file(3:4));



%% サンプリング周波数
[wdata,Fs] = audioread([path file]) ;
dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]
ymin=0; ymax=1.5;

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
T=tem1((S2(No).Var10+S2(No).Var11)/2);

%% グラフデータ
f1=figure(1);
set(f1,'Position', [500 100 1600 800])

%% グラフの切り取り範囲
if 1 <= rem(str2double(path_name1(1:2)),8) && rem(str2double(path_name1(1:2)),8) <=4
    timew=350;
else
    timew=150;
end

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

Amp=zeros(2,3);
Tp=zeros(2,3);


for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(pos+8).xpos,S1(pos+8).ypos,timew);
        
        figure(1);subplot(1,2,2)
        hold on
        plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','black')
        hold off
        xlim([4.4 4.55]);ylim([0 1.4]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title({strcat(num2str(No),path_name1(3:end));strcat(num2str(No+1),path_name1(3:end))});set(tp,tx)        
    end
    
    %% 反対側   
    file(3:4)=num2str(pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(oppos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(oppos+(j-1)*16,'%02d');
        file(9:10) = num2str(oppos+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(pos+8).xpos,S1(pos+8).ypos,timew);
        
        figure(1);subplot(1,2,2)
        hold on
        plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','black')
        hold off
    end
    
    %% 比較対象
    file(3:4) = num2str(17-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(9-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(9-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(9-pos+(j-1)*16,'%02d');
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(pos+8).xpos,S1(pos+8).ypos,timew);
        
        figure(1);subplot(1,2,2)
        hold on
        plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','blue')
        hold off

    end
    
    
    file(3:4)=num2str(9-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(17-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(17-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(17-pos+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(pos+8).xpos,S1(pos+8).ypos,timew);
        
        figure(1);subplot(1,2,2)
        hold on
        plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','red')
        hold off
    end  
end



%% from templeture to speed of sound
function T=tem1(t)
    T=331.5+0.6*t;
end



function[tp,amp]=denpan(iftr,d,delay,Fc,T,xtime,path,file,x1,y1,x2,y2,timew)
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
        idxmbase=fix((sqrt((x1-x2)^2+(y1-y2)^2)*192000/T+0.18*192));
        if idxmbase < 0
            idxmbase = 1;
        end
        wdmpoint=idxmbase+96;
        
        % 振幅最大値
        [tymax]=max(abs(wdm(idxmbase:idxmbase+timew)));
  
        % 伝搬時間
        tppos=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        tp=xtime(fix(tppos-0.18*192));

        
        amp=0;
        for k=0:9
            amp=amp+abs(wdm(wdmpoint+k));
        end
        amp=amp/10;
       
end
