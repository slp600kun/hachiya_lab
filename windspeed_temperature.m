%
% wave ファイルのプロット
%
%%
clearvars
close all

% スプレッドシート
filename1 = '6_sonar_tr4x4_pos_b41.csv';
T1 = readtable(filename1);
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


%% サンプリング周波数
[wdata,Fs] = audioread([path file]) ;
dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]



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

%% グラフの切り取り範囲
if 1 <= rem(str2double(path_name1(1:2)),8) && rem(str2double(path_name1(1:2)),8) <=4
    timew=350;
else
    timew=150;
end

%% グラフデータ
f1=figure(1);
set(f1,'Position', [500 400 1500 500])



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

w=zeros(2,1);
wlist=zeros(8,6);
v=zeros(8,3);
c=zeros(8,3);
tem=zeros(8,3);
micpos=[12;11;10;9;16;15;14;13];

for i=1:3
    for ix=1:8
        rec=micpos(ix);
        
        path(end-2:end-1)= num2str(ix+(i-1)*16,'%02d');
        path_name2(8:9)= num2str(ix+(i-1)*16,'%02d');
        file(3:4) = num2str(rec,'%02d');
        file(9:10) = num2str(ix+(i-1)*16,'%02d');
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
        idxmbase=fix((leng(S1(ix).xpos,S1(ix).ypos,S1(rec).xpos,S1(rec).ypos)*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end
        time.sppos(ix,i)=idxmbase;
        time.sp(ix,i)=xtime(idxmbase);
        time.ep(ix,i)=xtime(idxmbase+timew);

        % 振幅最大値
        [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+timew)));
        tpos=idxmbase+tpos;
        tmax.tppos(ix,i)=tpos;
        tmax.tpy(ix,i)=tymax;
        tmax.tp(ix,i)=xtime(tpos);

        % 伝搬時間
        pt.tppos(ix,i)=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        pt.tp(ix,i)=xtime(pt.tppos(ix,i))-0.18;

       %% 反対側
        path(end-2:end-1) = num2str(rec+(i-1)*16,'%02d');
        path_name2(8:9)= num2str(rec+(i-1)*16,'%02d');
        file(3:4) = num2str(ix,'%02d');
        file(9:10) = num2str(rec+(i-1)*16,'%02d');

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
        idxmbase=fix((leng(S1(ix).xpos,S1(ix).ypos,S1(rec).xpos,S1(rec).ypos)*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end

        time.smpos(ix,i)=idxmbase;
        time.sm(ix,i)=xtime(idxmbase);
        time.em(ix,i)=xtime(idxmbase+timew);

        % 振幅最大値
        [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+timew)));
        tpos=idxmbase+tpos;
        tmax.tmpos(ix,i)=tpos;
        tmax.tmy(ix,i)=tymax;
        tmax.tm(ix,i)=xtime(tpos);

        % 伝搬時間
        pt.tmpos(ix,i)=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        pt.tm(ix,i)=xtime(pt.tmpos(ix,i))-0.18;
    end
    
    %% 風速
    for iv=1:8
        v(iv,i)=wind(leng(S1(iv).xpos,S1(iv).ypos,S1(micpos(iv)).xpos,S1(micpos(iv)).ypos),pt.tp(iv,i)/1000,pt.tm(iv,i)/1000);
    end

    %% 風速ベクトル
    figure(1);
    subplot(1,3,i);
    
    hold on
    for iw=1:8
        ver=[S1(micpos(iw)).xpos-S1(iw).xpos ; S1(micpos(iw)).ypos-S1(iw).ypos];
        w=normc(ver).*v(iw,i);  %%風速を各方向の単位ベクトルに掛ける
        quiver((S1(iw).xpos+S1(micpos(iw)).xpos)/2,(S1(iw).ypos+S1(micpos(iw)).ypos)/2,w(1),w(2));
        axis ([-0.75 0.75 -0.75 0.75]);
        wlist(iw,2*i-1:2*i)=w(:);
    end    

    hold off
    
    %% 音速
    for ic=1:8
        c(ic,i)=sound(leng(S1(ic).xpos,S1(ic).ypos,S1(micpos(ic)).xpos,S1(micpos(ic)).ypos),pt.tp(ic,i),pt.tm(ic,i));
    end
    
    %% 温度
    for it=1:8
        tem(it,i)=tem2(1000*c(it,i));
    end    


end

disp(tem)

%% from templeture to speed of sound
function T=tem1(t)
    T=331.5+0.6*t;
end

%% from templeture to speed of sound
function T=tem2(t)
    T=(abs(t)-331.5)/0.6;
end

%% length
function len=leng(x1,y1,x2,y2)
    l=(x1-x2)^2+(y1-y2)^2;
    len=sqrt(l);
end

%% 風速
function v=wind(l,tp,tm)
    v=2*l*(tm-tp)/(tm+tp)^2;
end

%% 音速
function c=sound(l,tp,tm)
    c=2*l/(tp+tm);
end