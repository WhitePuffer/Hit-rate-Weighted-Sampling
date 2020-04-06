function [ ds,P] = distPointAffspace_All( L, X )
%DISTPOINTAFFSPACE compute the distance between a point X and an affine
%subspace L
     m = max(size(X));
     ds = zeros(m,1);
     for j=1:m
         ds(j) =  distPointAffspace_All_data( X(:,j) , L );
     end
     P = L(:);         
end

 function  [ d ] = distPointAffspace_All_data( X, L ) 
    %%% for one data toward on mode hypothesis
    f=size(X,1);
    L=reshape(L,f,4);
    % proietto  X su U ottenendo Y
    U=L(:,2:end); %giacitura di L
    p=L(:,1); % punto di L
    k=(X-p)'*U; %coefficienti di Y in U
    Y=p;
    for i=1:length(k)
         Y = Y + k(i).*U(:,i); %proiezione di X di L
    end
    d= norm(X-Y);
 end


