%%5-fold cross-validation


A = APDpart1(:,1:numPart*numData);
Ai = APDpart1(:,(numPart*numData+1):numPart*(numData+numDataI) );

% 创建5-fold cross-validation 分区
cv = cvpartition(size(A, 2), 'KFold', 5);

% 初始化存储训练集和测试集的元胞数组
TrainingB = cell(cv.NumTestSets, 1); % pure signal 训练集
TestC = cell(cv.NumTestSets, 1); % pure signal测试集

TrainingBi = cell(cv.NumTestSets, 1); % noise signal 训练集
TestCi = cell(cv.NumTestSets, 1); % pure signal 测试集

TrainingSet = cell(cv.NumTestSets, 1);
TestSet = cell(cv.NumTestSets, 1);
% 分割数据
for i = 1:cv.NumTestSets
    % 获取训练和测试索引
    trainIdx = training(cv, i); % 训练集索引
    testIdx = test(cv, i);      % 测试集索引
    
    % 分别提取训练集和测试集
    TrainingB{i} = A(:, trainIdx); % 训练集矩阵 (按列提取)
    TestC{i} = A(:, testIdx);  % 测试集矩阵 (按列提取)

    TrainingBi{i} = Ai(:, trainIdx); % 训练集矩阵 (按列提取)
    TestCi{i} = Ai(:, testIdx);  % 测试集矩阵 (按列提取)
end

% 显示每折训练集和测试集大小
for i = 1:cv.NumTestSets
    fprintf('Fold %d: Training Set Size = [%d, %d], Test Set Size = [%d, %d]\n', ...
        i, size(TrainingB{i}, 1), size(TrainingB{i}, 2), size(TestC{i}, 1), size(TestC{i}, 2));

    fprintf('Fold %d: Training Set Size = [%d, %d], Test Set Size = [%d, %d]\n', ...
        i, size(TrainingBi{i}, 1), size(TrainingBi{i}, 2), size(TestCi{i}, 1), size(TestCi{i}, 2));    
end
for i = 1:cv.NumTestSets
 TrainingSet{i}= [TrainingB{i} TrainingBi{i}];
 TestSet{i}= [TestC{i} TestCi{i}];
end
disp('3rd: End of Training Set B, Test Set C');
%TrainingB{1},TestC{1},2-5组成的。
%%  1-96列是正常信号
%%  96+1 - 96+96是含噪声信号