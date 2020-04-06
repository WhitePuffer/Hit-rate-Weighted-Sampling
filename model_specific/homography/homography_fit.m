% 单应性矩阵估计-直接线性变换(DLT)算法
% Homography fitting function using Direct Linear Transformation (DLT)
% 输入: 对应的两组齐次坐标 x1 x2; 每组 3 * N个元素 (x,y,w)'-[1~n]
% 输出: Homograph H 3*3
function H = homography_fit(X) 
%     [~,n] = size(A); % the size of matrix A  
%     [U,D,V] = svd(A);
%     x = V*(D(1:n,1:n)\U(:,1:n)'*b);
    x1 = X(1:3,:);
    x2 = X(4:6,:);
    [~,n] = size(x1); % n_x1 = n_x2
    A = zeros(2*n,9); % A is a 2n X 9 matrix
    Zr = [0 0 0]; 
    for i = 1:n
        A(2*i-1,:) = [Zr, -x2(3,i)*x1(:,i)', x2(2,i)*x1(:,i)'];
        A(2*i,:) = [x2(3,i)*x1(:,i)', Zr, -x2(1,i)*x1(:,i)'];
    end
    [Un,Dnn,V] = svd(A,0);  
    h = V(:,end);
    H = reshape(h,3,3)';
    H = H / H(3,3); % set H(3,3) -> 1 (Non-essential)
end
