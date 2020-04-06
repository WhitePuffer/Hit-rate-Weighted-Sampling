function [ fitfn, resfn, degenfn, psize, numpar,fit_model] = getModelParaXmu(model_type)

%---------------------------
% Model specific parameters.
%---------------------------
%-------------------------
%%% fitfn ֱ������������
%%% fit_model �����ͱ�ǩһ�����
%-------------------------
switch model_type 
    
    case 'line'
        fitfn = @line_fit;
        resfn = @line_res;
        degenfn = @line_degen;
        psize = 2;
        numpar = 3;       
        fit_model = @fit_lines_label;
    case 'circle'
        fitfn = @fittingfn_circle2d;
        resfn = @distfn_circle2d;
        degenfn = @degenfn_circle2d;
        psize = 3;
        numpar = 3;      
        fit_model = @fit_circles_label;
    case 'homography'
        fitfn = @homography_fit;
        resfn = @homography_res;
        degenfn = @homography_degen;
        psize = 6 %4+2;  %  6 
        numpar = 9;        
        fit_model = @fit_Homography_label;
    case 'fundamental'
        fitfn = @fundamental_fit;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 7;
        numpar = 9;      
    case 'fundamental8' %%% 8 �㷨
        fitfn = @fundamental_fit8;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 8;
        numpar = 9;      
        fit_model = @fit_Fundamental_label;
    case 'affspace_All'
        fitfn=@fit_aff;
        resfn=@distPointAffspace_All; 
        degenfn= @dummydegenerate;% @isAff3degenerate;
        psize = 6; %6;    %%% 4+2
        numpar = 24;  %%% ��������
        fit_model = nan;
    case 'MotionFourSubspace_155'
        fitfn=@fitplanebase;
        resfn=@planeptdist;
        degenfn=@isdegenerate;
        psize= 4; 
        numpar= 25; 
    case 'Subspace_Tip'
        fitfn = @motion_fit;
        resfn = @motion_res;
        degenfn = @motion_degen;
        psize = 6;
        numpar = 25; 
        fit_model = nan;
    case 'plane'
        fitfn=@plane_fit;
        resfn=@plane_res;
        degenfn=@plane_degen;        
        psize=3;
        numpar=4;       
        fit_model = @fit_Planes_label;       
    case 'affspace3'
        fitfn=@fit_aff;   %%%% fit_aff_label һ����
        resfn=@distPointAffspace; 
        degenfn= @dummydegenerate;% @isAff3degenerate;
        psize = 6;    %%% 4+2
        numpar = 6;
        fit_model=@fit_aff_label;
   case 'motion'
        fitfn = @motion_fit;
        resfn = @motion_res;
        degenfn = @motion_degen;
        psize = 4;
        numpar = 25;
        
       % d=nan;   %%% d  ��������        
       % hpFun=@hpAffspace;
 %--------- �����Ĳв�ͼ������ɵĸ�ʽ��ͬ �����ǩlabel �����---------------- 
   case 'segment'
        distFun = @distPointSeg;  %%% ��Ҫ�ٿ���ô��
        hpFun = @hpSeg;
        fit_model = @fit_seg;   %%% û�к�������
        psize = 2; 
        numpar =4;       
        isdegen= @dummydegenerate;
  case 'subspace4'
        distFun = @distPointSubspace4;
        hpFun = @hpSubspace4;
        fit_model = @fit_Subspace4;
        psize = 4;
        isdegen= @dummydegenerate;%@isSub4degenerate;
        numpar = nan;
  case 'subspace9'
        distFun = @distPointSubspace9;
        hpFun = @hpSubspace9;
        fit_model = @fit_Subspace9_label;
        cardmss = 9;
        isdegen= @dummydegenerate;
        d = nan;
 %------------ �����Ĳв�ͼ������ɵĸ�ʽ��ͬ ----------------
        
    otherwise
        error('unknown model type!');
end

end