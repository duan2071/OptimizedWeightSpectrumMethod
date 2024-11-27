%内容接u8_1.m,u8_2.m  分析w>0的部分，然后再分析APP
close all;
%% 3-----the Fourier transform of data
Fs=2e7;%20MHz
numD=numData;  % Number of pure signals
numI=numDataI; % Number of signals with interference

i_fold = 1;   %1-5
%% Spectral analysis of APD received signals
APDfreH=zeros(detd5/2+1,numPart*(numD+numI));
APDfreA=zeros(detd5,numPart*(numD+numI));
APDfreAll=zeros(detd5,numPart*(numD+numI));
for j=1:1:numPart*(numD+numI)*0.8
    APDy=TrainingSet{i_fold}(:,j);
    APDN=length(APDy);
    APDX=fft(APDy);       %FFT transform, the result exists in Y
    APDfreA(2:end,j)=abs(APDX(2:end) )/sum(abs(APDX(2:end)));   % Normalization, removing DC
    APDfreAll(:,j)=APDX;
    % Frequency scale of Y
    df=Fs/APDN;        % Frequency interval df
    tf = Fs*(0:(APDN/2))/APDN;  % frequency scale
    tfall= Fs*(0:(APDN-1))/APDN;  %All frequency scale
    % Correct the amplitude of the positive frequency
    APDX=APDX/APDN;
    z=abs(APDX);    % Amplitude of X

    % Extract the positive frequencies in X and merge the negative frequencies in X into the positive ones.
    APDY=APDX(1:APDN/2+1);               % Extract the positive frequency part of X.
    APDY(2:end-1) = 2*APDY(2:end-1);  % The portion of X with negative frequencies is merged into the portion with positive frequencies.

    % Calculate the amplitude and phase angle of the sequence Y in the frequency domain
    APDA=abs(APDY);       % Calculate the amplitude of Y

    pha=angle(APDY);   % Calculating the phase angle, The units are radians.

    APDA(1,1)=0;
    APDfreH(:,j)=APDA./sum(APDA);
end

%% Fig 4.
figure(101);
subplot(311);  plot((1:1:2048)*1e6/2e7,TrainingSet{i_fold}(:,60),'r');
ylim([1.1 1.7]);xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Normal signals');

subplot(312);  plot((1:1:2048)*1e6/2e7,TrainingSet{i_fold}(:,(96+60) ),'r');
ylim([1.2 1.5]);
xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Signals with interference');

figure(110);
subplot(311);  plot(tfall(1,2:end),real(APDfreA(2:end,60)),'r');
xticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]*1e7);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'});
xlabel('Frequency (MHz)');ylabel('Amplitude');title('Spectrum of normal signals');
subplot(312);  plot(tfall(1,2:end),real(APDfreA(2:end,(96+60) )),'r');
xticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]*1e7);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'});
xlabel('Frequency (MHz)');ylabel('Amplitude');title('Spectrum of signals with interference');

disp('第三步 傅里叶变换完成');


%% 4-----Data weighting analysis.
% parameter
    nP=numPart*numD*0.8;       % Number of pure signals
    nQ=numPart*numI*0.8;       % Number of signals with interference
    numNFS=zeros(detd5,(nP+nQ)); %2048 X 14
    
%   ****************The following 2 lines need to be modified*********************
    numNFS(:,1:nP)=APDfreA(:,1:nP);     %1-np is the pure signals, NF is the normalized data, Each 1 column is a set of data.
    numNFS(:,(nP+1):(nP+nQ) )=APDfreA(:,(nP+1):(nP+nQ) ); %1-10是含有噪声信号，放数据的前面。（接收噪声信号APD的频谱）
    numx1=ones(1,nP+nQ);     %1 rpw and N columns
    nX=[numNFS.' numx1.'].';   %transpose  % The last line is label.  1 column is a set of data.
 
    nN=length(numNFS(:,1));         %The number of frequency components of a set of data, i.e., the row vector of the NF.一组数据的频率成分数量，即NF的行向量
    %初始化起始点
    nbeta = 0*ones(nN+1,1);  %迭代变量的初始值设为0.   ----输入变量---参数：可修改
    nH=zeros(nP+nQ,1);
    for i=1:1:nP+nQ
        nH(i,1)=exp(dot((nbeta.'),nX(:,i)))./(1+exp(dot((nbeta.'),nX(:,i))) );   %待求  notes:dot点积函数
    end
    
    nY=[ones(1,nP) zeros(1,nQ)].';  %标签----正常数据标签=1，噪声数据标签=0
    nlambda=5;   %超参 衰减系数 如何设置影响衰减3

% 定义目标函数的梯度或损失函数 需要已知P,Q,X,H,Y,lambda
ngrad_f = @(beta) (1/2/(nP+nQ))*(nX*(nH-nY))+1/(nP+nQ)*nlambda.*beta;  %求偏导的公式

% Setup Parameters
nlearning_rate = 0.01; %Learning rate - parameter: modifiable.学习率---参数：可修改
nmax_iterations=100000;    %Setting the maximum number of iterations - parameter: modifiable. 设置最大迭代次数---参数：可修改
ntolerance=1e-6;     %Tolerance - Parameter: Modifiable.   容许误差---参数：可修改


%gradient descent formula
tic
for titeration = 1:nmax_iterations
    %Calculating the gradient   计算梯度
    ngradient =ngrad_f(nbeta);
    
    %Update parameters 更新参数
    nbeta = nbeta -nlearning_rate * ngradient;
    
    %Checking for convergence  检查收敛性
    if norm(ngradient) < ntolerance
        break;
    end
end
toc

% Show results 显示结果
fprintf('Optimal solution: beta = [ %f ]\n',nbeta);
fprintf('Number of iterations: %d \n',titeration);

% Counting the number of beta elements less than 0.  以下是统计beta元素小于0的个数。
nindex=find(nbeta(1:end-1,1)<=0);  %beta contains w b.    是w b的组合 只判断小于等于0的
num0 = numel(nindex); % Display number less than 0.       显示小于0的个数。
% 显示多少个OWS大于0，多少个OWS小于0.
% 本来写OWS，现在改为WS。
fprintf('WS <=0 Number:[ %f ]\n',num0);
fprintf('WS >0 Number:[ %f ]\n',length(nbeta(1:end-1,1))-num0);%
fprintf('WS >0 percentage:[ %f ]%\n',(length(nbeta(1:end-1,1))-num0)*100/length(nbeta));
wsP=length(nbeta(1:end-1,1))-num0;%权重大于0的数量
[~, ws_ordered] = sort(nbeta(1:end-1,1), 'descend');   %--将序号按照数值升序排列
wsser = ws_ordered(1:wsP,1);%保留的频率成分；
%生成一个序列
wsod = zeros(2048,1);
wsod(wsser,1)=1;
wsod(1,1)=1;  %WS保留的频率成分


% 根据结果，对比情况
% figure;
% plot(abs(nbeta(1:end-1,1)),'.-r');hold on;    %权重谱
% plot(abs(tzfreA(1:end,1)),'b');                %第1列数据，正常信号的频谱
% title('OWS vs 调制信号');legend('OWS','调制信号');
% 
% figure; % 第一个直流成分不看，从第二个开始看
% plot(jfall(1,2:end),nbeta(2:end-1,1),'r');hold on;
% plot(jfall(1,2:end),abs(APDfreA(2:end,1)),'b');
% title('OWS vs APD信号');legend('OWS','APD信号');

%% Fig 5.
figure(111);
subplot(211);  plot(tfall(1,1:end),nbeta(1:end-1,1),'m');
ylim([-0.02 0.12]);
xticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]*1e7);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'});
xlabel('Frequency (MHz)');ylabel('Amplitude');title('Optimized weight spectrum');

%% 5------整理出权重谱存在tzw中，tzwt中只保留正的部分。
nzw=nbeta(1:end-1,1);
nzwt=(nzw).'; %wt是1行n列
% 小于0的部分置零，并使用图像显示查看一下。
% figure;%对应图35
% subplot(2,1,1); plot(nzwt(2:end));hold on;  xlabel('Frequency');ylabel('w');title('raw w')

%nzwt(nzwt<0)=0; %beta的值小于0，则置为0
nzwt(1,nindex)=0; %小于0的成分置0
% subplot(2,1,2); plot(nzwt(2:end));   xlabel('Frequency');ylabel('w');title('The positive part of w')



%% 5.5--------------拟合再去除一部分
nn0=sum(nzwt~=0);     % wt不为0的个数，表示的是不为0的个数。
[~, wt_ordered] = sort(nzwt, 'descend');   %--将nzwt序号按照数值降序排列,wt_ordered存放是nzwt中排序的索引。
[~, ws_ordered2] = sort(nzwt, 'descend'); 
app = ones(nn0,2);%初始化向量
for i=1:nn0
    
    app(i,1)= nzwt(1,wt_ordered(1,i));         % 法向量的幅值排序，第一列存放法向量的值（或称权重值）,app第1列 存放了从大到小的幅值。
    app(i,2)= (nn0-i+1)/nn0;     % 计算百分比，第二列存放幅度值的百分比，分子是从小到大的排序，
%    app(i,2)=  app(i,1)/max(nzwt);%按最大值比例不好。
%        app(i,2)=  app(i,1)/sum(nzwt);%分母是大于0的幅度和，分子是幅度值，归一化处理。
end
    %app 1行是一个坐标（幅值，百分比），幅值是按大到小排列的。
figure %y轴是w的幅度，x仅仅表示序列，其对应的频率没有体现，对应的频率在app第2列中。
plot(app(:,1),app(:,2),'.'); hold on
title('Amplitude sequence');xlabel('sequence');ylabel('Amplitude of w');


%% step2-2 选取超参 ，即确定频谱中需要保留的幅值
% 拟合曲线的数量
xlen = length(find(nbeta>0));    %取幅值大于0的**************************
% 定义变量，存放均方差
SE1=zeros(xlen,1);        %保存第1段 误差平方的临时变量
RMSE1=zeros(xlen,1);    %第1段线性方程的误差平方和的开方
SE2=zeros(xlen,1);      %保存第2段 误差平方的临时变量
RMSE2=zeros(xlen,1);    %第2段线性方程的误差平方和的开方
RMSESUM=zeros(xlen,1);  %1+2段总体误差和
% 以下循环遍历拟合 
xs=1:1:xlen;xs=xs.';
for i=2:xlen-2
    c=polyfit(app(1:i,1),app(1:i,2),4); %c返回的是多项式的系数
    y=polyval(c,app(1:i,1));             %序号和幅值拟合多项式，y是拟合的结果
    ytest=app(1:i,2);       %已知的结果
    SE1(1:i,1)=abs(ytest-y).^2;    %求误差的平方
    RMSE1(i,1)=sqrt(sum(SE1(1:i,1)));%误差平方和，再开平方
    
    c2=polyfit(app(i+1:xlen,1),app(i+1:xlen,2),1); %c2返回的是多项式的系数
    y2=polyval(c2,app(i+1:xlen,1));     %序号和幅值拟合多项式，y是拟合的结果
    SE2(i+1:xlen,1)=abs(app(i+1:xlen,2)-y2).^2;%求误差的平方
    RMSE2(i,1)=sqrt(sum(SE2(i+1:xlen,1)));%误差平方和，再开平方
    RMSESUM(i,1)=RMSE1(i,1)+RMSE2(i,1);
end
% 拟合误差按降序排列
[~, min_ordered] = sort(RMSESUM, 'ascend');   %--将序号按照数值升序排列
% 误差最小的序号为k，1表示第一个不为0的数
k=find(RMSESUM(min_ordered(1:end)),1);% ***** 找到了最小幅值在数组app中的序号，第一个误差不为0的序列；k表示min_ordered的序列
anum=min_ordered(k);            %保留频率成分的数量84
% anum=166;   % 手动保留复制的个数========================
aRMSE=RMSESUM(min_ordered(k),1);%对应的误差
aAmp=app(anum,1);     %对应的幅值

% 找到的幅值 
fprintf(' Number of reserve:[ %f ]\n',anum);
%% 
% input para: nbeta，其中包含w和b的数量
% output: anum保留的频率成分的个数
%         aAmp是保留的最小的幅值（与anum对应）

%% Fig 7. fit a curve to a model.
figure(104);
subplot(211);plot(app(:,1),app(:,2),'ko', 'MarkerSize',4,'MarkerFaceColor','k');hold on;%点
%拟合的两个曲线
    nzc=polyfit(app(1:84,1),app(1:84,2),4); %c返回的是多项式的系数
    nzy=polyval(nzc,app(1:84,1));             %序号和幅值拟合多项式，y是拟合的结果
subplot(211);plot(app(1:84,1), nzy,'m','LineWidth',1.5);hold on;
    nzc2=polyfit(app(84+1:xlen,1),app(84+1:xlen,2),1); %c2返回的是多项式的系数
    nzy2=polyval(nzc2,app(84+1:xlen,1));     %序号和幅值拟合多项式，y是拟合的结果
subplot(211);plot(app(85:end,1), nzy2,'r','LineWidth',1.5);hold on;
xlabel('Normal vector');ylabel('Percentage of sequences');

% close all;
%% 6------去除不要频谱成分的幅值
nwt3noise=nzwt(1:end-1)';%权重谱保留的部分, nwzt is the part of w>0.
[~, ws_ordered_ows] = sort(nwt3noise, 'descend');%按降序排列
%% OWS保留的频率成分
owsod =zeros(2048,1);
owsod(ws_ordered_ows(1:anum,1),1)=1;
owsod(1,1)=1;% ows保留的成分

% %6.1-----APP拟合保留幅度值的内容
for i=1:1:length(nwt3noise)
    if(nwt3noise(i)<aAmp)
        nwt3noise(i)=0; % Zero for small amplitude
    end
end
nsigf=find(nwt3noise==0); %   保存不要频率成分的序号(w=0对应的序号)


%% 6------去除不要频谱成分的幅值
nIXnoise2 = zeros(length(jdfreAll(:,1)),nP+nQ); %每一列是预测的数据，行是2048个数据
for j=1:1:nP+nQ
    % 
    % Spectrum of received signals.    
    tf8ton=APDfreAll(:,j); % 2N X 1 matrix. 2048X1

    tf7ton=tf8ton;

    tf7ton(nsigf(2:end,1),1) = 0; % The excess frequency components are cleared to zero, but the first DC component is retained.

tf3ton=tf7ton;


    %% ifft
    % in wt2noise is the frequency component of the reservation. wt2noise中是保留的频率成分
    nIXnoise2(:,j)=ifft(tf3ton(:,1)); % the Fourier inverse transform.

end

%% Fig 6
figure(102);
subplot(321); plot((1:1:2048)*1e6/2e7,tzpart1(:,60),'b');hold on;       %xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('transmitted signal (1)');axis tight
subplot(323);  plot((1:1:2048)*1e6/2e7, TrainingB{i_fold}(:,60),'r');hold on;
ylim([1.1 1.7]);    %xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Normal signals');axis tight
subplot(325)        %
plot((1:1:2048)*1e6/2e7,real(nIXnoise2(:,60)),'m');hold on;
ylim([1.1 1.7]);    %xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Normal signals separated by OWS');axis tight


subplot(322); plot((1:1:2048)*1e6/2e7,tzpart1(:,198),'b');hold on;xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('transmitted signal (2)');axis tight

subplot(324);  plot((1:1:2048)*1e6/2e7,TrainingB{i_fold}(:,198),'r');
%ylim([1.2 1.5]);
%xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Signals with interference');axis tight
subplot(326) %
plot((1:1:2048)*1e6/2e7,real(nIXnoise2(:,198)),'m');hold on;
ylim([1.3 1.5]);
%xlim([0 110]);
xlabel('Time(\mu s)');ylabel('Amplitude');title('Signals with interference separated by OWS');axis tight


disp('u8_4 Run Completion');
