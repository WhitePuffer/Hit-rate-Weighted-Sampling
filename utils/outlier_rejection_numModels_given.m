% modified by jianglan
function [ C_new ] = outlier_rejection_numModels_given( C,  numModels)
%OUTREJHIST outlier rejection based on cardinality: the rationale is that
%outliers tends to emerge as small cluster in conceptual space.
%T-Linkage is agonostic about the outlier rejection strategy adopted: some
%improvement can be obtained exploiting the randomness of a model (code will be available).
% INPUT:
%       C: clustering vector
%       cardmss: size of minimal sample set
%OUTPUT
%       C_new: clustering vector, outlier are labeled 0.
%
%
% Author: Luca Magri 2014
% For any comments, questions or suggestions about the code please contact
% luca (dot) magri (at) unimi (dot) it


    label=sort(unique(C));
    N=length(label); 
    card=zeros(1,N);
    for l=1:N
        card(l)=sum(C==label(l));
    end

    label=[label;];
    [~, v]=sort(card,'ascend'); %升序
    label_sorted=label(v);

    cut = N-numModels;
    C_new = C;
    for j=1:cut
        C_new(C_new==label_sorted(j))=0;
    end
    [C_new(C_new~=0),~,~] = grp2idx(C_new(C_new~=0));
end

