function [ M ] =fit_aff_label(X,C)
%fit_aff fit a 3 dimensional affine space


cardmss=4;
f=size(X,1);

if  ~exist('C','var') || isempty(C)
    
    [m,n]= size(X);
    centroid=sum(X,2)./n;
    % momentum matrix
    momentum_matrix=zeros(m);
    for j=1:n
        momentum_matrix=momentum_matrix+(X(:,j)-centroid)*(X(:,j)-centroid)';
    end
    [U,~,~]=svd(momentum_matrix);
    U=U(:,1:cardmss-1);
    U=U(:);
    
    M=vertcat(centroid, U);
    
    
else
    
    
    
    label=unique(C);
    N=length(label); %numero di cluster;
    
    M=zeros(f*cardmss,N);
    
    for i=1:N
        
        L  = label(i);
        points2fit = X(:,C==L);
        [m,n]= size(points2fit);
        
        %centroid;
        centroid=sum(points2fit,2)./n;
        % momentum matrix
        momentum_matrix=zeros(m);
        
        for j=1:n
            momentum_matrix=momentum_matrix+(points2fit(:,j)-centroid)*(points2fit(:,j)-centroid)';
        end
        [U,~,~]=svd(momentum_matrix);
        U=U(:,1:cardmss-1);
        U=U(:);
        %
        M(:,i)=vertcat(centroid, U);
    end
    
end
