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
micpos = [12;11;10;9;16;15;14;13];


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
set(f1,'Position', [500 100 1000 500])
%% グラフデータ
f2=figure(2);
set(f2,'Position', [600 100 1000 500])

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

%% 凡例
tiledlayout(2,1)
str1=' 2 \rightarrow 11 ';
str2=' 3 \rightarrow 10 ';
str3=' 11 \rightarrow 2 ';
str4=' 10 \rightarrow 3 ';
str5=' 6 \rightarrow 15 ';
str6=' 7 \rightarrow 14 ';
str7=' 15 \rightarrow 6 ';
str8=' 14 \rightarrow 7 ';
strA={str1,str6,str3,str8};
strB={str2,str5,str4,str7};
title1="通常時の経路間波形伝搬時間と振幅(ヒートガン直上でない)";
title2="通常時の経路間波形伝搬時間と振幅(ヒートガン直上)";
title3="ヒートガン設置時の経路間波形伝搬時間と振幅(ヒートガン直上でない)";
title4="ヒートガン設置時の経路間波形伝搬時間と振幅(ヒートガン直上)";


for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    %% 2往路
    file(3:4)=num2str(micpos(pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a1(i,j)=Amp(i,j);
        amp.stdA(1)=std(amp.a1,1,'all');
        
        figure(2);subplot(1,2,1)
        hold on
        p1 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','black'); set(p1,ps)
        xlim([4.35 4.5]);ylim([0 1]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title(title3);set(tp,tx)        
    end
    %% 7往路
    file(3:4) = num2str(micpos(9-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(9-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(9-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(9-pos+(j-1)*16,'%02d');
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a2(i,j)=Amp(i,j);
        amp.stdA(2)=std(amp.a2,1,'all');
        figure(2);subplot(1,2,1)
        p6 =plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','magenta'); set(p6,ps)

    end        

    
    %% 2復路 
    file(3:4)=num2str(pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a3(i,j)=Amp(i,j);
        amp.stdA(3)=std(amp.a3,1,'all');
        figure(2);subplot(1,2,1)
        p3 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','black'); set(p3,ps)
    end   
    %% 7復路
    file(3:4)=num2str(9-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(9-pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a4(i,j)=Amp(i,j);
        amp.stdA(4)=std(amp.a4,1,'all');
        figure(2);subplot(1,2,1)
        p8 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','magenta'); set(p8,ps)
    end        
   %% 3往路
   file(3:4)=num2str(micpos(pos+1),'%02d');
   for j=1:3
        path(end-2:end-1)= num2str(pos+1+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+1+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+1+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a5(i,j)=Amp(i,j);
        amp.stdA(5)=std(amp.a5,1,'all');
        figure(2);subplot(1,2,2)
        hold on
        p2 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','red'); set(p2,ps)
        xlim([4.35 4.5]);ylim([0 1]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title(title4);set(tp,tx)
    end    
    %% 6往路
    file(3:4) = num2str(micpos(8-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(8-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(8-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(8-pos+(j-1)*16,'%02d');
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a6(i,j)=Amp(i,j);
        amp.stdA(6)=std(amp.a6,1,'all');
        figure(2);subplot(1,2,2)
        p5 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','blue'); set(p5,ps)
        tp=title(title4);set(tp,tx)
    end

    %% 3復路
    file(3:4)=num2str(pos+1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos+1)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a7(i,j)=Amp(i,j);
        amp.stdA(7)=std(amp.a7,1,'all');
        figure(2);subplot(1,2,2)
        p4 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','red'); set(p4,ps)
    end     
    %% 6復路
    file(3:4)=num2str(8-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(8-pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.a8(i,j)=Amp(i,j);
        amp.stdA(8)=std(amp.a8,1,'all');
        figure(2);subplot(1,2,2)
        p7 =plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','blue'); set(p7,ps)
    end
 
end
figure(1);subplot(1,2,2)
legend([p1 p6 p3 p8],strA,'Location','northeast');


figure(2);subplot(1,2,2)
legend([p2 p5 p4 p7],strB,'Location','northeast');

%% デフォルト
path = erase(path,path(pos_wn+4:pos_frame-2));
subpos_wn = strfind(path_name1,'wn');
path_name1 = erase(path_name1,path_name1(subpos_wn+4:end));
for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(i,'%02d');
    path_name1(1:2) = num2str(i,'%02d');
    file(3:4)=num2str(micpos(pos),'%02d');
    %% 2往路
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b1(i,j)=Amp(i,j);
        amp.stdB(1)=std(amp.b1,1,'all');
        
        figure(1);subplot(1,2,1)
        hold on
        p1 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','black'); set(p1,ps)
        xlim([4.35 4.5]);ylim([0 1]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title(title1);set(tp,tx)    
    end
    %% 7往路
    file(3:4) = num2str(micpos(9-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(9-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(9-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(9-pos+(j-1)*16,'%02d');
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b2(i,j)=Amp(i,j);
        amp.stdB(2)=std(amp.b2,1,'all');
        figure(1);subplot(1,2,1)
        p6 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','magenta'); set(p6,ps)

    end     
    %% 2復路 
    file(3:4)=num2str(pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b3(i,j)=Amp(i,j);
        amp.stdB(3)=std(amp.b3,1,'all');
        figure(1);subplot(1,2,1)
        p3 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','black'); set(p3,ps)
    end
    %% 7復路
    file(3:4)=num2str(9-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(9-pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b4(i,j)=Amp(i,j);
        amp.stdB(4)=std(amp.b4,1,'all');
        figure(1);subplot(1,2,1)
        p8 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','magenta'); set(p8,ps)
    end      
   %% 3往路
   file(3:4)=num2str(micpos(pos+1),'%02d');
   for j=1:3
        path(end-2:end-1)= num2str(pos+1+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+1+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+1+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b5(i,j)=Amp(i,j);
        amp.stdB(5)=std(amp.b5,1,'all');
        figure(1);subplot(1,2,2)
        hold on
        p2 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','red'); set(p2,ps)     
        xlim([4.35 4.5]);ylim([0 1]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title(title2);set(tp,tx)        
   end    
    %% 6往路
    file(3:4) = num2str(micpos(8-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(8-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(8-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(8-pos+(j-1)*16,'%02d');
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b6(i,j)=Amp(i,j);
        amp.stdB(6)=std(amp.b6,1,'all');
        figure(1);subplot(1,2,2)
        p5 = plot(Tp(i,j),Amp(i,j),'o','MarkerEdgeColor','blue'); set(p5,ps)
    end
    %% 3復路
    file(3:4)=num2str(pos+1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos+1)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b7(i,j)=Amp(i,j);
        amp.stdB(7)=std(amp.b7,1,'all');
        figure(1);subplot(1,2,2)
        p4 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','red'); set(p4,ps)
    end      
    %% 6復路
    file(3:4)=num2str(8-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(8-pos)+(j-1)*16,'%02d');
        
        [Tp(i,j),Amp(i,j)]=denpan(iftr,d,delay,Fc,T,xtime,path,file,S1(pos).xpos,S1(pos).ypos,S1(micpos(pos)).xpos,S1(micpos(pos)).ypos,timew);
        amp.b8(i,j)=Amp(i,j);
        amp.stdB(8)=std(amp.b8,1,'all');
        figure(1);subplot(1,2,2)
        p7 = plot(Tp(i,j),Amp(i,j),'s','MarkerEdgeColor','blue'); set(p7,ps)
    end
   
end
figure(1);subplot(1,2,1)
legend([p1 p6 p3 p8],strA,'Location','northeast');


figure(2);subplot(1,2,1)
legend([p2 p5 p4 p7],strB,'Location','northeast');




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
        idxmbase=fix((sqrt((x1-x2)^2+(y1-y2)^2)*192000/T-50+0.18*192));
        if idxmbase < 0
            idxmbase = 1;
        end
        wdmpoint=fix(idxmbase+146);
        
        % 振幅最大値
        [tymax]=max(abs(wdm(idxmbase:idxmbase+timew)));
  
        % 伝搬時間
        tppos=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.05), 1 )+idxmbase;
        tp=xtime(tppos)-0.18;

        
        amp=0;
        for k=0:9
            amp=amp+abs(wdm(+wdmpoint+k));
        end
        amp=amp/10;
       
       
end
