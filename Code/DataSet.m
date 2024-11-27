%% 按列划分
close all;clc;clear;

% 创建50列数据矩阵A (随机数据示例, 50行x10列)
rng(42); % 固定随机种子
A = rand(1024, 50);

% 创建5-fold cross-validation 分区
cv = cvpartition(size(A, 2), 'KFold', 5);

% 初始化存储训练集和测试集的元胞数组
B = cell(cv.NumTestSets, 1); % 训练集
C = cell(cv.NumTestSets, 1); % 测试集

% 分割数据
for i = 1:cv.NumTestSets
    % 获取训练和测试索引
    trainIdx = training(cv, i); % 训练集索引
    testIdx = test(cv, i);      % 测试集索引
    
    % 分别提取训练集和测试集
    B{i} = A(:, trainIdx); % 训练集矩阵 (按列提取)
    C{i} = A(:, testIdx);  % 测试集矩阵 (按列提取)
end

% 显示每折训练集和测试集大小
for i = 1:cv.NumTestSets
    fprintf('Fold %d: Training Set Size = [%d, %d], Test Set Size = [%d, %d]\n', ...
        i, size(B{i}, 1), size(B{i}, 2), size(C{i}, 1), size(C{i}, 2));
end

