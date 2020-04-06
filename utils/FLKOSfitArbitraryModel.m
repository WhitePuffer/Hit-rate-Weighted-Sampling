function [theta_f, sigma_f, Cinl, ht,  f_inx] = FLKOSfitArbitraryModel(x, k, model_type,Threshold,W,filename)
% INPUTS
%x = data - D x N matrix of data; D - dimention of data ; N - Number of data points
%k - The size of the minimum acceptible structure in application
%model_type - type of model to fit
    %'line2D', 'plane3D', 'homography', 'fundamental', 'subspace'
%Threshold = For dicotomizing inliers from outliers: 
    %T = 2~3.0 assuming normally distribured noise;

% OUTPUTS
%theta_f - final set of paramers after one FLKOS iteration 模型参数f
%sigma_f - final sigma 尺度
%Cinl - indexes of the inlier to theta_f 内点index 
%ht - the residuals with respect to theta_f 点到f的残差
%f_inx - the sample that produced theta_f 生成f的采样集


[ fitfn, resfn, degenfn, psize, numpar ] = getModelParam_Pplus2(model_type); %获取模型参数
[~,N] = size(x); % number of points
n_rand_inits = 1; %number of random initializations 初始化最多次数
n_iterations = 50; %number of flkos iterations 迭代最大次数

BestJ = 1e7; %初始化最佳值 cost function

init_counter = 0;
while(init_counter < n_rand_inits)
    
    %initlialize from a random psize-tupple 初始: 随机抽取p的点 p:=p+2 (比最小采样集个数大2)
    [psub,pinx] = datasample(x,psize,2, 'Replace', false, 'Weights', W);
%     plabel = zeros(1,N);
%     plabel(pinx) = 1;
%     subplot(2,1,1);
%     plot_img_demo(filename, plabel);

%     pinx = randsample(n,psize);
%     psub = x(:,pinx);
    iter_counter = 0; % 迭代次数计数
    while( (iter_counter < n_iterations)  )
        theta = feval(fitfn,psub);  % Fit the model on the p-subset 拟合
        dist = feval(resfn,theta,x);% Compute residuals. 计算残差
        
        [SRes,I] = sort(dist); %残差排序
        % 获取LkOS , 残差排序为(k-psize-1):k 共psize个点
        pinx = I( (k-psize+1):k ); % Get new tupple around kth sorted residual 
        psub = x(:,pinx);
        
        % 计算LkOS的残差和
        Jnew = sum( SRes( (k-psize+1):k ) ); % cost func - sum of p residuals(squared) around K
        
        if (Jnew < BestJ) % store the best value 如果(残差和)比当前最佳好,则更新当前最佳
            BestJ = Jnew;
            theta_f = theta;
            f_inx = pinx;
        end
        if (iter_counter ==0)
            old_pinx1 = pinx;
            old_pinx = pinx;
        end
        % old_pinx1 为第i-2次的最佳LkOS ; old_pinx 为第i-1次的最佳LkOS ;
        % dist:停止条件计算
        dist = feval(resfn,theta,x(:,old_pinx)); 
        dist = mean(dist);
        dist1 = feval(resfn,theta,x(:,old_pinx1));
        dist1 = mean(dist1);
        
        sd = SRes(k);% the kth residual with respect to current theta 当前第k大残差
        
        %stopping criterion
        %check if the tupples in last two interations are still within the 
        %k residuals
        if ( ( dist < sd ) && ( dist1 < sd )  &&  (iter_counter > 10)) %停止条件判断
%            %这里有加了一个trick, 必须超过10次迭代才判断停止条件
             break;
        end
        old_pinx1 = old_pinx;
        old_pinx = pinx;
   
        iter_counter = iter_counter+1;
    end
%     if (Jnew < BestJ) % store the best value
%             BestJ = Jnew;
%             theta_f = theta;
%             f_inx = pinx;
%         end
    init_counter = init_counter+1;
end

%get the best theta and refine it 最佳假设refine
ht = feval(resfn,theta_f,x);
[sresd2,sindx]=sort(ht);
sigma_f = 0.01;
Cinl = sresd2< 2 * sigma_f;
% [Cinl,sigma_f]  = getInliers(sresd2,sindx,k,Threshold);
% [ sigma_f, Cinl ] = Fast_AVG_MSSE( ht, Threshold, 100, 0, 2 );
%get the inliers based  on refined theta
% psub = x(:,Cinl);  
% theta_f = feval(fitfn,psub); 
% ht = feval(resfn,theta_f,x);
% [sresd2,sindx]=sort(ht);
% [Cinl, sigma_f] = getInliers(sresd2,sindx,k,Threshold);
% psub = x(:,Cinl);  
% theta_f = feval(fitfn,psub); 
end
