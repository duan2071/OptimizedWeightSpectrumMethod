%���ݽ�u8_1.m,u8_2.m  ����w>0�Ĳ��֣�Ȼ���ٷ���APP
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

disp('������ ����Ҷ�任���');


%% 4-----Data weighting analysis.
% parameter
    nP=numPart*numD*0.8;       % Number of pure signals
    nQ=numPart*numI*0.8;       % Number of signals with interference
    numNFS=zeros(detd5,(nP+nQ)); %2048 X 14
    
%   ****************The following 2 lines need to be modified*********************
    numNFS(:,1:nP)=APDfreA(:,1:nP);     %1-np is the pure signals, NF is the normalized data, Each 1 column is a set of data.
    numNFS(:,(nP+1):(nP+nQ) )=APDfreA(:,(nP+1):(nP+nQ) ); %1-10�Ǻ��������źţ������ݵ�ǰ�档�����������ź�APD��Ƶ�ף�
    numx1=ones(1,nP+nQ);     %1 rpw and N columns
    nX=[numNFS.' numx1.'].';   %transpose  % The last line is label.  1 column is a set of data.
 
    nN=length(numNFS(:,1));         %The number of frequency components of a set of data, i.e., the row vector of the NF.һ�����ݵ�Ƶ�ʳɷ���������NF��������
    %��ʼ����ʼ��
    nbeta = 0*ones(nN+1,1);  %���������ĳ�ʼֵ��Ϊ0.   ----�������---���������޸�
    nH=zeros(nP+nQ,1);
    for i=1:1:nP+nQ
        nH(i,1)=exp(dot((nbeta.'),nX(:,i)))./(1+exp(dot((nbeta.'),nX(:,i))) );   %����  notes:dot�������
    end
    
    nY=[ones(1,nP) zeros(1,nQ)].';  %��ǩ----�������ݱ�ǩ=1���������ݱ�ǩ=0
    nlambda=5;   %���� ˥��ϵ�� �������Ӱ��˥��3

% ����Ŀ�꺯�����ݶȻ���ʧ���� ��Ҫ��֪P,Q,X,H,Y,lambda
ngrad_f = @(beta) (1/2/(nP+nQ))*(nX*(nH-nY))+1/(nP+nQ)*nlambda.*beta;  %��ƫ���Ĺ�ʽ

% Setup Parameters
nlearning_rate = 0.01; %Learning rate - parameter: modifiable.ѧϰ��---���������޸�
nmax_iterations=100000;    %Setting the maximum number of iterations - parameter: modifiable. ��������������---���������޸�
ntolerance=1e-6;     %Tolerance - Parameter: Modifiable.   �������---���������޸�


%gradient descent formula
tic
for titeration = 1:nmax_iterations
    %Calculating the gradient   �����ݶ�
    ngradient =ngrad_f(nbeta);
    
    %Update parameters ���²���
    nbeta = nbeta -nlearning_rate * ngradient;
    
    %Checking for convergence  ���������
    if norm(ngradient) < ntolerance
        break;
    end
end
toc

% Show results ��ʾ���
fprintf('Optimal solution: beta = [ %f ]\n',nbeta);
fprintf('Number of iterations: %d \n',titeration);

% Counting the number of beta elements less than 0.  ������ͳ��betaԪ��С��0�ĸ�����
nindex=find(nbeta(1:end-1,1)<=0);  %beta contains w b.    ��w b����� ֻ�ж�С�ڵ���0��
num0 = numel(nindex); % Display number less than 0.       ��ʾС��0�ĸ�����
% ��ʾ���ٸ�OWS����0�����ٸ�OWSС��0.
% ����дOWS�����ڸ�ΪWS��
fprintf('WS <=0 Number:[ %f ]\n',num0);
fprintf('WS >0 Number:[ %f ]\n',length(nbeta(1:end-1,1))-num0);%
fprintf('WS >0 percentage:[ %f ]%\n',(length(nbeta(1:end-1,1))-num0)*100/length(nbeta));
wsP=length(nbeta(1:end-1,1))-num0;%Ȩ�ش���0������
[~, ws_ordered] = sort(nbeta(1:end-1,1), 'descend');   %--����Ű�����ֵ��������
wsser = ws_ordered(1:wsP,1);%������Ƶ�ʳɷ֣�
%����һ������
wsod = zeros(2048,1);
wsod(wsser,1)=1;
wsod(1,1)=1;  %WS������Ƶ�ʳɷ�


% ���ݽ�����Ա����
% figure;
% plot(abs(nbeta(1:end-1,1)),'.-r');hold on;    %Ȩ����
% plot(abs(tzfreA(1:end,1)),'b');                %��1�����ݣ������źŵ�Ƶ��
% title('OWS vs �����ź�');legend('OWS','�����ź�');
% 
% figure; % ��һ��ֱ���ɷֲ������ӵڶ�����ʼ��
% plot(jfall(1,2:end),nbeta(2:end-1,1),'r');hold on;
% plot(jfall(1,2:end),abs(APDfreA(2:end,1)),'b');
% title('OWS vs APD�ź�');legend('OWS','APD�ź�');

%% Fig 5.
figure(111);
subplot(211);  plot(tfall(1,1:end),nbeta(1:end-1,1),'m');
ylim([-0.02 0.12]);
xticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]*1e7);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'});
xlabel('Frequency (MHz)');ylabel('Amplitude');title('Optimized weight spectrum');

%% 5------�����Ȩ���״���tzw�У�tzwt��ֻ�������Ĳ��֡�
nzw=nbeta(1:end-1,1);
nzwt=(nzw).'; %wt��1��n��
% С��0�Ĳ������㣬��ʹ��ͼ����ʾ�鿴һ�¡�
% figure;%��Ӧͼ35
% subplot(2,1,1); plot(nzwt(2:end));hold on;  xlabel('Frequency');ylabel('w');title('raw w')

%nzwt(nzwt<0)=0; %beta��ֵС��0������Ϊ0
nzwt(1,nindex)=0; %С��0�ĳɷ���0
% subplot(2,1,2); plot(nzwt(2:end));   xlabel('Frequency');ylabel('w');title('The positive part of w')



%% 5.5--------------�����ȥ��һ����
nn0=sum(nzwt~=0);     % wt��Ϊ0�ĸ�������ʾ���ǲ�Ϊ0�ĸ�����
[~, wt_ordered] = sort(nzwt, 'descend');   %--��nzwt��Ű�����ֵ��������,wt_ordered�����nzwt�������������
[~, ws_ordered2] = sort(nzwt, 'descend'); 
app = ones(nn0,2);%��ʼ������
for i=1:nn0
    
    app(i,1)= nzwt(1,wt_ordered(1,i));         % �������ķ�ֵ���򣬵�һ�д�ŷ�������ֵ�����Ȩ��ֵ��,app��1�� ����˴Ӵ�С�ķ�ֵ��
    app(i,2)= (nn0-i+1)/nn0;     % ����ٷֱȣ��ڶ��д�ŷ���ֵ�İٷֱȣ������Ǵ�С���������
%    app(i,2)=  app(i,1)/max(nzwt);%�����ֵ�������á�
%        app(i,2)=  app(i,1)/sum(nzwt);%��ĸ�Ǵ���0�ķ��Ⱥͣ������Ƿ���ֵ����һ������
end
    %app 1����һ�����꣨��ֵ���ٷֱȣ�����ֵ�ǰ���С���еġ�
figure %y����w�ķ��ȣ�x������ʾ���У����Ӧ��Ƶ��û�����֣���Ӧ��Ƶ����app��2���С�
plot(app(:,1),app(:,2),'.'); hold on
title('Amplitude sequence');xlabel('sequence');ylabel('Amplitude of w');


%% step2-2 ѡȡ���� ����ȷ��Ƶ������Ҫ�����ķ�ֵ
% ������ߵ�����
xlen = length(find(nbeta>0));    %ȡ��ֵ����0��**************************
% �����������ž�����
SE1=zeros(xlen,1);        %�����1�� ���ƽ������ʱ����
RMSE1=zeros(xlen,1);    %��1�����Է��̵����ƽ���͵Ŀ���
SE2=zeros(xlen,1);      %�����2�� ���ƽ������ʱ����
RMSE2=zeros(xlen,1);    %��2�����Է��̵����ƽ���͵Ŀ���
RMSESUM=zeros(xlen,1);  %1+2����������
% ����ѭ��������� 
xs=1:1:xlen;xs=xs.';
for i=2:xlen-2
    c=polyfit(app(1:i,1),app(1:i,2),4); %c���ص��Ƕ���ʽ��ϵ��
    y=polyval(c,app(1:i,1));             %��źͷ�ֵ��϶���ʽ��y����ϵĽ��
    ytest=app(1:i,2);       %��֪�Ľ��
    SE1(1:i,1)=abs(ytest-y).^2;    %������ƽ��
    RMSE1(i,1)=sqrt(sum(SE1(1:i,1)));%���ƽ���ͣ��ٿ�ƽ��
    
    c2=polyfit(app(i+1:xlen,1),app(i+1:xlen,2),1); %c2���ص��Ƕ���ʽ��ϵ��
    y2=polyval(c2,app(i+1:xlen,1));     %��źͷ�ֵ��϶���ʽ��y����ϵĽ��
    SE2(i+1:xlen,1)=abs(app(i+1:xlen,2)-y2).^2;%������ƽ��
    RMSE2(i,1)=sqrt(sum(SE2(i+1:xlen,1)));%���ƽ���ͣ��ٿ�ƽ��
    RMSESUM(i,1)=RMSE1(i,1)+RMSE2(i,1);
end
% �������������
[~, min_ordered] = sort(RMSESUM, 'ascend');   %--����Ű�����ֵ��������
% �����С�����Ϊk��1��ʾ��һ����Ϊ0����
k=find(RMSESUM(min_ordered(1:end)),1);% ***** �ҵ�����С��ֵ������app�е���ţ���һ����Ϊ0�����У�k��ʾmin_ordered������
anum=min_ordered(k);            %����Ƶ�ʳɷֵ�����84
% anum=166;   % �ֶ��������Ƶĸ���========================
aRMSE=RMSESUM(min_ordered(k),1);%��Ӧ�����
aAmp=app(anum,1);     %��Ӧ�ķ�ֵ

% �ҵ��ķ�ֵ 
fprintf(' Number of reserve:[ %f ]\n',anum);
%% 
% input para: nbeta�����а���w��b������
% output: anum������Ƶ�ʳɷֵĸ���
%         aAmp�Ǳ�������С�ķ�ֵ����anum��Ӧ��

%% Fig 7. fit a curve to a model.
figure(104);
subplot(211);plot(app(:,1),app(:,2),'ko', 'MarkerSize',4,'MarkerFaceColor','k');hold on;%��
%��ϵ���������
    nzc=polyfit(app(1:84,1),app(1:84,2),4); %c���ص��Ƕ���ʽ��ϵ��
    nzy=polyval(nzc,app(1:84,1));             %��źͷ�ֵ��϶���ʽ��y����ϵĽ��
subplot(211);plot(app(1:84,1), nzy,'m','LineWidth',1.5);hold on;
    nzc2=polyfit(app(84+1:xlen,1),app(84+1:xlen,2),1); %c2���ص��Ƕ���ʽ��ϵ��
    nzy2=polyval(nzc2,app(84+1:xlen,1));     %��źͷ�ֵ��϶���ʽ��y����ϵĽ��
subplot(211);plot(app(85:end,1), nzy2,'r','LineWidth',1.5);hold on;
xlabel('Normal vector');ylabel('Percentage of sequences');

% close all;
%% 6------ȥ����ҪƵ�׳ɷֵķ�ֵ
nwt3noise=nzwt(1:end-1)';%Ȩ���ױ����Ĳ���, nwzt is the part of w>0.
[~, ws_ordered_ows] = sort(nwt3noise, 'descend');%����������
%% OWS������Ƶ�ʳɷ�
owsod =zeros(2048,1);
owsod(ws_ordered_ows(1:anum,1),1)=1;
owsod(1,1)=1;% ows�����ĳɷ�

% %6.1-----APP��ϱ�������ֵ������
for i=1:1:length(nwt3noise)
    if(nwt3noise(i)<aAmp)
        nwt3noise(i)=0; % Zero for small amplitude
    end
end
nsigf=find(nwt3noise==0); %   ���治ҪƵ�ʳɷֵ����(w=0��Ӧ�����)


%% 6------ȥ����ҪƵ�׳ɷֵķ�ֵ
nIXnoise2 = zeros(length(jdfreAll(:,1)),nP+nQ); %ÿһ����Ԥ������ݣ�����2048������
for j=1:1:nP+nQ
    % 
    % Spectrum of received signals.    
    tf8ton=APDfreAll(:,j); % 2N X 1 matrix. 2048X1

    tf7ton=tf8ton;

    tf7ton(nsigf(2:end,1),1) = 0; % The excess frequency components are cleared to zero, but the first DC component is retained.

tf3ton=tf7ton;


    %% ifft
    % in wt2noise is the frequency component of the reservation. wt2noise���Ǳ�����Ƶ�ʳɷ�
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
