function [ par,res,inx ] = ProximateSample( data, method, M, sigma, model_type )
% m is the number of hypothesis, cardmss is the cand of this model.
if(nargin < 2)
    return;
end

dataForDistance = data(1:2,:);

% if(nargin <= 5)
    K = pdist(dataForDistance', 'euclidean');%����������֮��ľ���,�洢��������K��
    D = squareform(K);%��������Kת���������ǻ��������Ǿ�����ʽ,�洢��D��
    K = [];
% end

switch(method)
    case 'max'
        sortedD = sort(D(1:size(dataForDistance,2),1:size(dataForDistance,2)),1,'ascend');
        range = max(sortedD(par1,:));
        nearPtsTab = zeros(size(D,1), size(D,2));
        nearPtsTab(find(D < range)) = 1;
        minpwmean = mean(sortedD(2,:));
    case 'min'
        sortedD = sort(D(1:size(dataForDistance,2),1:size(dataForDistance,2)),1,'ascend');
        range = min(sortedD(par1,:));
        nearPtsTab = zeros(size(D,1), size(D,2));
        nearPtsTab(find(D < range)) = 1;
        minpwmean = mean(sortedD(2,:));
    case 'mean'
        sortedD = sort(D(1:size(dataForDistance,2),1:size(dataForDistance,2)),1,'ascend');
        range = mean(sortedD(par1,:));
        nearPtsTab = zeros(size(D,1), size(D,2));
        nearPtsTab(find(D < range)) = 1;
        minpwmean = mean(sortedD(2,:));
    case 'maxkout'
       sortedD = sort(D(1:size(dataForDistance,2),1:size(dataForDistance,2)),1,'ascend');
       sortCol = sort(sortedD(par1,:),'descend');
        range = sortCol(par2);
        nearPtsTab = zeros(size(D,1), size(D,2));
        nearPtsTab(find(D < range)) = 1;
        minpwmean = mean(sortedD(2,:));
    case 'exp'
        nearPtsTab = exp(-(D.^2)/sigma);
    otherwise
        nearPtsTab = exp(-(D.^2)/sigma);
end

for i=1:size(D,1)
    nearPtsTab(i,i) = 0;
    nearPtsTab(i,:) = nearPtsTab(i,:) / sum(nearPtsTab(i,:));   
    nearPtsTab(i,:) = cumsum(nearPtsTab(i,:));   
    %��ÿһ�н��й�һ������������㵽�õ�ľ����һ��
end

[ fitfn,resfn,degenfn,psize,numpar ] = getModelParam(model_type);
%---------------------����M��ģ��--------------------------------
n = size(data,2);
par = zeros(numpar,M);
res = zeros(n,M);
inx = zeros(psize,M);
trialcount = 1;
maxDataTrials = 50;
%   hw  =   waitbar(0,'Generating Hypotesis ...');
while (trialcount <= M)
%     waitbar(trialcount/maxTrials,hw);
    % Select at random s datapoints to form a trial model, M.
    degenerate = 1;
    count = 1;
    while degenerate && count <= maxDataTrials;
      % Generate 1 random index in the range 1..npts
      ind = randsample(n, 1);
      nearPointsCdf = nearPtsTab(ind,:);
      rnd = rand(psize,1);
      [dum, ind] = histc(rnd,[0 nearPointsCdf]);
      if ind==0 
          continue;
      end
      curModSample = data(:,ind);
      degenerate = 1;
      if(length(ind) == length(unique(ind)))
        degenerate = feval(degenfn, curModSample);
      end
      % Test that these points are not a degenerate configuration.
      	    
      if ~degenerate;
        % Fit model to this random selection of data points.
        
        tempPar = feval(fitfn, curModSample);
      end

      count = count + 1;
    end

    if degenerate;
      warning ( 'MATLAB:ransac:Output', ...
                'Unable to select a nondegenerate data set!' );
      break;
    end

    res(:,trialcount) = feval(resfn, tempPar, data);
    tempPar = reshape(tempPar, numpar, 1);%%change the matrix shape 
    par(:,trialcount) = tempPar;
    inx(:,trialcount) = ind;
    trialcount = trialcount + 1;  
end
%-------------------------------------------------------------------
end


