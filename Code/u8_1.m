%求基波信号和调制信号
clc;clear;close all;
%% 1-----加载文件目录
%% APD信号目录
%1-10正常的APD信号，11-20是含噪的APD信号
filename_APD_ch1_1='..\Dataset\scope_0_1.csv';%1
filename_APD_ch1_2='..\Dataset\scope_1_1.csv';%2
filename_APD_ch1_3='..\Dataset\scope_2_1.csv';%3
filename_APD_ch1_4='..\Dataset\scope_3_1.csv';%4
filename_APD_ch1_5='..\Dataset\scope_4_1.csv';%5
filename_APD_ch1_6='..\Dataset\scope_5_1.csv';%1
filename_APD_ch1_7='..\Dataset\scope_6_1.csv';%7
filename_APD_ch1_8='..\Dataset\scope_7_1.csv';%8
filename_APD_ch1_9='..\Dataset\scope_8_1.csv';%9
filename_APD_ch1_10='..\Dataset\scope_9_1.csv';%10

filename_APD_ch1_11='..\Dataset\scope_10_1.csv';%11
filename_APD_ch1_12='..\Dataset\scope_11_1.csv';%12



%% 读出文件
APD_ch1_1 = readmatrix(filename_APD_ch1_1);%1
APD_ch1_2 = readmatrix(filename_APD_ch1_2);%2
APD_ch1_3 = readmatrix(filename_APD_ch1_3);%3
APD_ch1_4 = readmatrix(filename_APD_ch1_4);%4
APD_ch1_5 = readmatrix(filename_APD_ch1_5);%5
APD_ch1_6 = readmatrix(filename_APD_ch1_6);%6
APD_ch1_7 = readmatrix(filename_APD_ch1_7);%7
APD_ch1_8 = readmatrix(filename_APD_ch1_8);%8
APD_ch1_9 = readmatrix(filename_APD_ch1_9);%9
APD_ch1_10 = readmatrix(filename_APD_ch1_10);%10

APD_ch1_11 = readmatrix(filename_APD_ch1_11);%11
APD_ch1_12 = readmatrix(filename_APD_ch1_12);%12



%% 仅读出APD信号
%数据量numData
numData=6;
numDataI=6;
APD=zeros(2e6,numData+numDataI);

APD(:,1)= APD_ch1_1(3:length(APD_ch1_1(:,2)),2:2);% 
APD(:,2)= APD_ch1_2(3:length(APD_ch1_2(:,2)),2:2);% 
APD(:,3)= APD_ch1_3(3:length(APD_ch1_3(:,2)),2:2);% 
APD(:,4)= APD_ch1_4(3:length(APD_ch1_4(:,2)),2:2);% 
APD(:,5)= APD_ch1_5(3:length(APD_ch1_5(:,2)),2:2);% 
APD(:,6)= APD_ch1_6(3:length(APD_ch1_6(:,2)),2:2);% 
APD(:,7)= APD_ch1_7(3:length(APD_ch1_7(:,2)),2:2);% 
APD(:,8)= APD_ch1_8(3:length(APD_ch1_8(:,2)),2:2);% 
APD(:,9)= APD_ch1_9(3:length(APD_ch1_9(:,2)),2:2);% 
APD(:,10)= APD_ch1_10(3:length(APD_ch1_10(:,2)),2:2);%

APD(:,11)= APD_ch1_11(3:length(APD_ch1_11(:,2)),2:2);% 
APD(:,12)= APD_ch1_12(3:length(APD_ch1_12(:,2)),2:2);% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调制信号目录
filename_tz_ch2_1='..\Dataset\scope_0_2.csv';%1
filename_tz_ch2_2='..\Dataset\scope_1_2.csv';%2
filename_tz_ch2_3='..\Dataset\scope_2_2.csv';%3
filename_tz_ch2_4='..\Dataset\scope_3_2.csv';%4
filename_tz_ch2_5='..\Dataset\scope_4_2.csv';%5
filename_tz_ch2_6='..\Dataset\scope_5_2.csv';%1
filename_tz_ch2_7='..\Dataset\scope_6_2.csv';%7
filename_tz_ch2_8='..\Dataset\scope_7_2.csv';%8
filename_tz_ch2_9='..\Dataset\scope_8_2.csv';%9
filename_tz_ch2_10='..\Dataset\scope_9_2.csv';%10
% 以下是含有噪声的接收信号
filename_tz_ch2_11='..\Dataset\scope_10_2.csv';%11
filename_tz_ch2_12='..\Dataset\scope_11_2.csv';%12
filename_tz_ch2_13='..\Dataset\scope_12_2.csv';%13



%% 读出文件
tz_ch2_1 = readmatrix(filename_tz_ch2_1);%1
tz_ch2_2 = readmatrix(filename_tz_ch2_2);%2
tz_ch2_3 = readmatrix(filename_tz_ch2_3);%3
tz_ch2_4 = readmatrix(filename_tz_ch2_4);%4
tz_ch2_5 = readmatrix(filename_tz_ch2_5);%5
tz_ch2_6 = readmatrix(filename_tz_ch2_6);%6
tz_ch2_7 = readmatrix(filename_tz_ch2_7);%7
tz_ch2_8 = readmatrix(filename_tz_ch2_8);%8
tz_ch2_9 = readmatrix(filename_tz_ch2_9);%9
tz_ch2_10 = readmatrix(filename_tz_ch2_10);%10
tz_ch2_11 = readmatrix(filename_tz_ch2_11);%11
tz_ch2_12 = readmatrix(filename_tz_ch2_12);%12



%% 仅读出调制信号
tz=zeros(2e6,numData+numDataI);

tz(:,1)= tz_ch2_1(3:length(tz_ch2_1(:,2)),2:2);% 
tz(:,2)= tz_ch2_2(3:length(tz_ch2_2(:,2)),2:2);% 
tz(:,3)= tz_ch2_3(3:length(tz_ch2_3(:,2)),2:2);% 
tz(:,4)= tz_ch2_4(3:length(tz_ch2_4(:,2)),2:2);% 
tz(:,5)= tz_ch2_5(3:length(tz_ch2_5(:,2)),2:2);% 
tz(:,6)= tz_ch2_6(3:length(tz_ch2_6(:,2)),2:2);% 
tz(:,7)= tz_ch2_7(3:length(tz_ch2_7(:,2)),2:2);% 
tz(:,8)= tz_ch2_8(3:length(tz_ch2_8(:,2)),2:2);% 
tz(:,9)= tz_ch2_9(3:length(tz_ch2_9(:,2)),2:2);% 
tz(:,10)= tz_ch2_10(3:length(tz_ch2_10(:,2)),2:2);% 

tz(:,11)= tz_ch2_11(3:length(tz_ch2_11(:,2)),2:2);% 
tz(:,12)= tz_ch2_12(3:length(tz_ch2_12(:,2)),2:2);% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 基带信号目录
filename_jidai_ch3_1='..\Dataset\scope_0_3.csv';%1
filename_jidai_ch3_2='..\Dataset\scope_1_3.csv';%2
filename_jidai_ch3_3='..\Dataset\scope_2_3.csv';%3
filename_jidai_ch3_4='..\Dataset\scope_3_3.csv';%4
filename_jidai_ch3_5='..\Dataset\scope_4_3.csv';%5
filename_jidai_ch3_6='..\Dataset\scope_5_3.csv';%1
filename_jidai_ch3_7='..\Dataset\scope_6_3.csv';%7
filename_jidai_ch3_8='..\Dataset\scope_7_3.csv';%8
filename_jidai_ch3_9='..\Dataset\scope_8_3.csv';%9
filename_jidai_ch3_10='..\Dataset\scope_9_3.csv';%10

filename_jidai_ch3_11='..\Dataset\scope_10_3.csv';%11
filename_jidai_ch3_12='..\Dataset\scope_11_3.csv';%12



%% 读出基带信号对应文件
jidai_ch3_1 = readmatrix(filename_jidai_ch3_1);%1
jidai_ch3_2 = readmatrix(filename_jidai_ch3_2);%2
jidai_ch3_3 = readmatrix(filename_jidai_ch3_3);%3
jidai_ch3_4 = readmatrix(filename_jidai_ch3_4);%4
jidai_ch3_5 = readmatrix(filename_jidai_ch3_5);%5
jidai_ch3_6 = readmatrix(filename_jidai_ch3_6);%6
jidai_ch3_7 = readmatrix(filename_jidai_ch3_7);%7
jidai_ch3_8 = readmatrix(filename_jidai_ch3_8);%8
jidai_ch3_9 = readmatrix(filename_jidai_ch3_9);%9
jidai_ch3_10 = readmatrix(filename_jidai_ch3_10);%10

jidai_ch3_11 = readmatrix(filename_jidai_ch3_11);%11
jidai_ch3_12 = readmatrix(filename_jidai_ch3_12);%12


%% 仅读出基带信号
jidai=zeros(2e6,numData+numDataI);

jidai(:,1)= jidai_ch3_1(3:length(jidai_ch3_1(:,2)),2:2);% 
jidai(:,2)= jidai_ch3_2(3:length(jidai_ch3_2(:,2)),2:2);% 
jidai(:,3)= jidai_ch3_3(3:length(jidai_ch3_3(:,2)),2:2);% 
jidai(:,4)= jidai_ch3_4(3:length(jidai_ch3_4(:,2)),2:2);% 
jidai(:,5)= jidai_ch3_5(3:length(jidai_ch3_5(:,2)),2:2);% 
jidai(:,6)= jidai_ch3_6(3:length(jidai_ch3_6(:,2)),2:2);% 
jidai(:,7)= jidai_ch3_7(3:length(jidai_ch3_7(:,2)),2:2);% 
jidai(:,8)= jidai_ch3_8(3:length(jidai_ch3_8(:,2)),2:2);% 
jidai(:,9)= jidai_ch3_9(3:length(jidai_ch3_9(:,2)),2:2);% 
jidai(:,10)= jidai_ch3_10(3:length(jidai_ch3_10(:,2)),2:2);% 

jidai(:,11)= jidai_ch3_11(3:length(jidai_ch3_11(:,2)),2:2);% 
jidai(:,12)= jidai_ch3_12(3:length(jidai_ch3_12(:,2)),2:2);% 



%APD每行是1组APD接收到的信号，tz每行是1组发射的调制信号，jidai每行是需要发射的基带信号
disp('1st: End of read.');





