function [ error,newLabel ] = errorForTransMatrix( groundLabel,resultLabel )
%ERRORFORTRANSMATRIX Summary of this function goes here
%   Detailed explanation goes here
    uniqueGround = unique(groundLabel); % groudtruth 结构index (含outliers)
    uniqueResult = unique(resultLabel); % resultLabel 结构index (含outliers)
    newLabel = zeros(size(resultLabel)); % 重标签 resultLabel
    nGround = length(uniqueGround); %  groudtruth 结构个数 (含outliers)
    nResult = length(uniqueResult); % resultLabel 结构个数 (含outliers)
    similarity = zeros(nGround,nResult); % 相似度 g与r各类之间的
    correct = 0;
    correspondence = zeros(1,nResult); % g与r的对应
    for i = 1:nGround
        for j = 1:nResult
            similarity(i,j) = length(intersect(find(groundLabel==uniqueGround(i)),find(resultLabel==uniqueResult(j))));
        end
    end
    for i = 1:min(nGround,nResult)
        maxSimilarity = max(max(similarity));
        [x,y] = find(similarity == maxSimilarity);
        correct = correct + maxSimilarity;
        correspondence(y) = x;
        similarity(x,:) = 0;
        similarity(:,y) = 0 ;
    end
    error = 1-correct/length(groundLabel);
%     fprintf('the misclassification points:%d\n',length(groundLabel)-correct);
    if(nResult>nGround)
        count = nGround;
        for i = 1:nResult
            if(correspondence(i) == 0)
                count = count + 1;
                correspondence(i) = count;
            end
        end
    end
    for i = 1:nResult
        newLabel(find(resultLabel==uniqueResult(i))) = correspondence(i)-1;
    end
    
end

