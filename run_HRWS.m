%% Initalization
close all;
clear;clc;
%------------------ Parameters
Threshold = 3.5; %Fundamental Scale Threshold
% sigma_f = 0.03; %for fundamental
% sigma_f = 0.04;%62; % for homo
sigma_f = 0.1% for SRNU
num_repetitions = 1;%20%50; %每张图 重复100次 记录数据
numHyp = 1000; %number of hypothesis to be generated for clustering
% [ fitfn, resfn, degenfn, psize, numpar ] = getModelPara_Xmu('homography');
% [ fitfn, resfn, degenfn, psize, numpar ] = getModelPara_Xmu('fundamental8');
[ fitfn, resfn, degenfn, psize, numpar ] = getModelPara_Xmu('affspace_All');
h = (Threshold^2) * (sigma_f^2);
%------------------ Path
addpath(genpath('./model_specific'));
addpath('./utils');
addpath('./utils/GMM/');
data_path = 'data/SNU_Based_Dataset_XMU_P_19_9_1/';
% data_path = 'data/AdelaiderRMF/fundamental/';
% data_path = 'data/AdelaiderRMF/homography/';
addpath(data_path);
data_files = dir(data_path);
data_files(1:2) = [];
result_path = ['Results_Affi_' ...
    'hyp' num2str(numHyp) '_'...
    'rep' num2str(num_repetitions) '_'...
    'AWKSH/'];
if(~exist(result_path ,'dir'))
    mkdir(result_path)
end
if(~exist([result_path 'TotalResults'],'dir'))
    mkdir([result_path 'TotalResults'])
end

%% Dataset Evaluation
seg_errors = zeros(length(data_files), num_repetitions);
seg_times = zeros(length(data_files), num_repetitions);
diary([result_path 'TotalResults/log.txt']); diary on;
for f=3%:length(data_files)
    filename = data_files(f).name;
    fprintf('Processing data %s: %s \n', num2str(f), data_files(f).name);
    pack = load(data_files(f).name);
%     data = pack.data;
%     label = pack.label;
    data = pack.xy;
    label = pack.SNU_GT;
    % ---------remove repeating rows in data-------
    [data,ia,~] = unique(data','rows');
    data = data';
    N = size(data,2);
    label = label(ia);
    numModels=max(label) - min(label);
    % --------- parameters setting -------
    k=floor(12); % 设置LkOS中的k值
    SampFrac = min([1/2, 1/(numModels)]) ; % 初始采样Ns个样本占总点数的百分比 Ns = N/numModels
    %---------- normalise2dpts(Golden Standard Algorithm)
    dat_img_1 = normalise2dpts(data(1:3,:));
    dat_img_2 = normalise2dpts(data(4:6,:));
    X = [dat_img_1;dat_img_2];
    %----------Cost Based Sampling------------
    sample_hitrate = zeros(1, N);
    hypo_hitrate = zeros(1, N);
    csp = zeros(1,numModels); % Clear Sample Point set
    H = zeros(N,numHyp);
    T = zeros(1,50);
    for i = 1:num_repetitions
        csp = zeros(1,numModels); % Clear Sample Point set
        sample_hitrate = zeros(1, N);
        hypo_hitrate = zeros(1, N);
        W = ones(1,N);
        t1 = tic;
        CI = zeros(2*numHyp, N);
        NN = N;
        for j = 1:numHyp
            % ----------- sampling ----------
            %             [Xs,Is] = datasample(X, floor(SampFrac*N), 2, 'Replace', false, 'Weights', W);
            [Xs,Is] = datasample(X, min(floor(N/2),NN), 2, 'Replace', false, 'Weights', W);
            sample_hitrate(Is) = sample_hitrate(Is)+1;
            Ws = W(Is);
            % ----------- Sample Generation ------------
            [~,Ns] = size(Xs); % number of points
            n_rand_inits = 1; %number of random initializations 初始化最多次数
            n_iterations = 50; %number of flkos iterations 迭代最大次数 原始为50
            init_counter = 0;
            pinliers_old = 0;
            update_counter = 0;
            pinliers_Best = 0;
            while(init_counter < n_rand_inits)
                %initlialize from a random psize-tupple 初始: 随机抽取p的点 p:=p+2 (比最小采样集个数大2)
                iter_counter = 0; % 迭代次数计数
                JBest = 0;
                [psub,psub_index] = datasample(Xs, psize, 2, 'Replace', false, 'Weights', Ws);
                ww = zeros(1,50);
                while(iter_counter < n_iterations)
                    theta = feval(fitfn,psub);  % Fit the model on the p-subset 拟合
                    res = feval(resfn,theta,Xs);% Compute residuals. 计算残差
                    [sorted_res, sorted_res_index] = sort(res); %残差排序
                    pinliers = sorted_res < h;
                    if sum(pinliers) >= k
                        Jnew = sum(1 - (sorted_res(pinliers)/h).^2)/sum(pinliers);
                        ww(iter_counter+1) = Jnew;
                        if (Jnew > JBest)
                            JBest = Jnew;
                            theta_f = theta;
                            update_counter = iter_counter;
                            f_inx = psub_index;
%                             if iter_counter >25
%                             end
                        end
                    end
                    if iter_counter >= 15 && iter_counter - update_counter >= 10
%                         a = sum(pinliers_old);
%                         b = sum(pinliers);
%                         c = sum(pinliers_old & pinliers);
%                         T(iter_counter+1) = T(iter_counter+1)+1;
                        break;
                    end
                    psub_index = sorted_res_index( (k-psize+1):k ); % Get new tupple around kth sorted residual
                    psub = Xs(:,psub_index); % 获取下一步迭代模型的LkOS , 残差排序为(k-psize-1):k 共psize个点
                    pinliers_old = pinliers;
                    iter_counter = iter_counter+1;
                end
                init_counter = init_counter+1;
            end                 
            label_check = label(Is(f_inx));
            if all(label_check == label_check(1) & label_check(1) ~= 0)
               csp(label_check(1)) = csp(label_check(1))+1;
            end
            % Sub-Sampling refine
            %---Best Hyp---
            res = feval(resfn,theta_f,Xs);
            theta_f = feval(fitfn,Xs(:,res < h));
            % Weight Adjustment 
            ht = feval(resfn, theta_f, X);
            Cinl = ht< h; % inliers 为 1, outliers 为 0 的逻辑数组
            Coutl = ~Cinl;
            hypo_hitrate(Cinl) = hypo_hitrate(Cinl)+1;
            CI(j, :) = Cinl;
            H(:,j) = exp(-ht/( 2*(sigma_f^2)));
            
            
            %---Newest Hyp----
%             res = feval(resfn,theta,Xs);
%             theta = feval(fitfn,Xs(:,res < h));
%             % Weight Adjustment 
%             ht2 = feval(resfn, theta, X);
%             Cinl2 = ht2< h; % inliers 为 1, outliers 为 0 的逻辑数组
%             Coutl2 = ~Cinl2;
%             hypo_hitrate(Cinl2) = hypo_hitrate(Cinl2)+1;
%             CI(2*j, :) = Cinl2;
%             H(:,2*j) = exp(-ht2/( 2*(sigma_f^2)));
            
            % -------------- update the bootstraping weights
%             W(W == 0) = 1;
            if (j == floor(numHyp/4))
                [keep, ~, th, gmm] = GMMremove(hypo_hitrate/j, 0.5);
                th2 = norminv(0.75, gmm.mu(2), gmm.sig(2));
            end
            if (j >= numHyp/4)
                W(hypo_hitrate  < (th*j) ) = 0.0001;
                
                W(hypo_hitrate  > (th2*j) ) =0.0001;
                NN = sum(hypo_hitrate >= th*j);
            end
            if numModels ~= 1
                W(Cinl) = W(Cinl)/2; % w<-w x 2; w(C_inl)<-w(C_inl)/4
                W(Coutl) = W(Coutl)*2;
%                 W(Cinl2) = W(Cinl2)/2; % w<-w x 2; w(C_inl)<-w(C_inl)/4
%                 W(Coutl2) = W(Coutl2)*2;
                W(Is) = W(Is)/2;
                W(~Is) = W(~Is)*2;
%                 W(W<.1) = .1; % 防止权重过小, [trick: 论文中没有]
                W(W>10) = 10; % w(w>w*) <- 1/N ; 其中w* = 20/N;
            else
                W = ones(1,N);
            end
        end

        % ---------- spectral Clustering ---------
%         G = H*H'; % 亲和矩阵G = H*H^T 生成
%         %         numModels = numModels-1;
%         [~, ClustLabels] = spectralClustering_ALI(G, numModels); %利用G 谱聚类ALI
%         % ---------- Evaluation ---------
%         time = toc(t1);
%         seg_times(f,i) = time;
%         [error,index] = missclass(ClustLabels-1,label);
%         seg_errors(f, i) = error;
%         disp([num2str(i) ' | ' 'time: ' num2str(time) ' | ' 'error: ' num2str(error)]);
%         disp(csp)
%         plot_img_demo(filename, ClustLabels);
%         pause; 
    end
%     fprintf('Mean error = %f \n', mean(seg_errors(f,:)));
%     % Sigle result save
%     Avgtime = mean(seg_times(f,:));
%     meanError = mean(seg_errors(f,:));
%     stdError = std(seg_errors(f,:));
%     medianError = median(seg_errors(f,:));
%     minError = min(seg_errors(f,:));
%     maxError = max(seg_errors(f,:));
%     save([result_path data_files(f).name], 'Avgtime', 'meanError', ...
%         'stdError', 'medianError', 'minError', 'maxError');
end
diary off;
%% Total Results Save
% [mean_err, Avg_times] = ResultSave(seg_times, seg_errors, result_path);
% fprintf('Total Mean error = %f \n', mean(mean_err));
% fprintf('Total Mean Time = %f \n', mean(Avg_times));
a = sample_hitrate/numHyp/num_repetitions;
b = hypo_hitrate/numHyp/num_repetitions;
hitrate_plot(a, b, label)
hhr = hypo_hitrate/numHyp/num_repetitions;