function [inliers,sigma] = getInliers_sigma(rs,I,k,sigma)
    inliers = I(rs< 2.5 * sigma);
end
    
