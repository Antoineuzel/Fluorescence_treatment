function [S385total,S405total,S385_corrected,S405_corrected,fluorophore,lambda,res_385_correc_Kim_full_exp,res_405_correc_Kim_full_exp,residu] = corrected_fluo_LS_WL_1(path,file,name_spectralon,isSmallSpectralon,min_wl_reflectance,max_wl_reflectance,min_wl_fluo,max_wl_fluo,fluorophores_385,fluorophores_405)
    %%%%% function to correct fluo when the WL acquisition is performed on
    %%%%% less wavelength compare to the fluorescence acquisition, typically when
    %%%%% acq is performed with LS WL 1, written by Antoine Uzel 24/08/2024
    n_medium=1.4;
    file_spectralon_theo='D:\Lyon thèse\data\theoretical data\spectralon\Reflectance_values_array.txt';

    %%Process calculation
    scale=load("D:\Lyon thèse\soft\setupcalibration\scale_new_theory_intralipids.mat");
    file_WL=strrep(file,'fluo','diffuseReflectances');
    [WLShort,WLLong,lambda_wl,WLight,idx_min,idx_max]=load_reflectance(fullfile(path,file_WL),min_wl_reflectance,max_wl_reflectance);
    if ~isfield(WLight,'lambda') 
    WLight.lambda=WLight.raw_lambda;
    end
    scale_short=interp1(scale.Scale.Short.WL,scale.Scale.Short.Sy,lambda_wl);
    n_short=interp1(scale.Scale.Short.WL,scale.Scale.Short.n,WLight.lambda);
     if isfield(WLight,'name_spectralon')
        file_spectralon_exp=fullfile(path,WLight.name_spectralon);
    else 
        file_spectralon_exp=fullfile(path,name_spectralon);
     end
    
    [SsourceShort,SsourceLong]=calibration_spectralon(file_spectralon_exp,file_spectralon_theo);
    
    
    if isSmallSpectralon % Use of small spectralon
        SsourceShort=SsourceShort/0.8;
    end
        R_Short=WLShort./(SsourceShort(idx_min:idx_max).*scale_short);
  % if one wants to fit reflectance  
    data_blood=importdata('D:\Lyon thèse\data\theoretical data\HbO2 Hb visible.txt');
    data_blood=data_blood.data;
    fat=load('D:\Lyon thèse\soft\soft manip Arthur\read txt files\mua_fat.mat');
    mua_fat=fat.absorbeur.mua;
    lambda_fat=fat.absorbeur.lambda;
    data_cyt_b_oxi=importdata("D:\Lyon thèse\data\theoretical data\cyt b oxidized.txt");
    data_cyt_b_red=importdata("D:\Lyon thèse\data\theoretical data\cyt b reduced.txt");
lambda_blood=data_blood(:,1);
lambda_cit_b_oxi=data_cyt_b_oxi(:,1);
cit_b_oxi=data_cyt_b_oxi(:,2);
lambda_cit_b_red=data_cyt_b_red(:,1);
cit_b_red=data_cyt_b_red(:,2);

HbO2=data_blood(:,2);
Hb=data_blood(:,3);

% Definition of blood constants useful later
CHb=150 ; %150 g/L théorique 
MHb=64500; % masse molaire de l'hémoglobine en g/mol
ua_Hb=2.303*Hb/MHb;
ua_HbO2=2.303*HbO2/MHb;
% Definition of the laser excitation wavelengths 
HbO2_interp=interp1(lambda_blood,ua_HbO2,WLight.lambda);
Hb_interp=interp1(lambda_blood,ua_Hb,WLight.lambda);
cit_b_oxi_interp=interp1(lambda_cit_b_oxi,cit_b_oxi,WLight.lambda);
cit_b_red_interp=interp1(lambda_cit_b_red,cit_b_red,WLight.lambda);
    % Extraction of scattering & absorption properties
lb = [0,0,0,0,0,0,0];
ub = [100,4,1,70,1,100];
    options=optimoptions('lsqcurvefit',...
        'MaxFunctionEvaluations',1e10,...
        'FunctionTolerance',1e-10, ...
        'MaxIter',1e8, ...
        'StepTolerance',1e-10);
    x0=[10,1,0,0.9,0.5,0.5,1,1];
    rho_short=240e-4;% in meters
    ua=[HbO2_interp;Hb_interp;mua_fat';cit_b_oxi_interp;cit_b_red_interp];

    fun_short=@(x,wl)reflectance_Kim(x,rho_short,ua,wl,n_medium,n_short,min_wl_reflectance,max_wl_reflectance);
    [x_short,residu.reflectance]= lsqcurvefit(fun_short,x0,WLight.lambda,R_Short,lb,ub,options);
    %% obtain fluo data
    [S385total,S405total,lambda]=load_fluo(fullfile(path,file),min_wl_fluo,max_wl_fluo);
%     R_short_extracted=reflectance_Kim(x_short,rho_short,ua,WLight.lambda,n_medium,n_short,min_wl_fluo,max_wl_fluo);
    R_short_extracted=R_Short;
    R_short_extracted=1;
    lb = [0,0,0,0,0,0,0];
    ub = [1,1000,1000,1000,1000,1000,1000];
    
    options=optimoptions('lsqcurvefit',...
        'MaxFunctionEvaluations',1e10,...
        'FunctionTolerance',1e-6, ...
        'MaxIter',1e8, ...
        'StepTolerance',1e-10); 

    
    
    fit_405=@(x,lambda)fluo_exp_385(x,fluorophores_405,lambda);
    fit_385=@(x,lambda)fluo_exp_385(x,fluorophores_385,lambda);
    lb = [0,0,0,0,0,0,0,494,14,589,9,636,5.5,618,7.5,0];
    ub = [1,1000,1000,100,1000,1000,1000,496,16,591,11,638,7.5,620,9,1000];

    x0=[0.1,0.1,0.1,0.01,0.1,0.1,0.1,495,10,590,10,619,6,636,6];
    [res_385_correc_Kim_full_exp,residu.fit_385] = lsqcurvefit(@(a,x) fit_385(a(2:end),x).*(R_short_extracted.^a(1)),x0,lambda,S385total,lb,ub,options);
    [res_405_correc_Kim_full_exp,residu.fit_405] = lsqcurvefit(@(a,x) fit_405(a(2:end),x).*(R_short_extracted.^a(1)),x0,lambda,S405total,lb,ub,options);
    res_385=res_385_correc_Kim_full_exp(2:end);
    res_405=res_405_correc_Kim_full_exp(2:end);

    S385_corrected=S385total./(R_short_extracted.^res_385_correc_Kim_full_exp(1));
    S405_corrected=S405total./(R_short_extracted.^res_405_correc_Kim_full_exp(1));
     
    fluorophore.flavine_385_exp=fit_385([res_385(1),0,0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.flavine_405_exp=fit_405([res_405(1),0,0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.NADH_385_exp=fit_385([0,res_385(2),0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.NADH_405_exp=fit_405([0,res_405(2),0,0,0,0,1,1,1,1,1,1,1,1],lambda);
    fluorophore.gaussian_385=fit_385([0,0,res_385(3),0,0,0,res_385(7),res_385(8),1,1,1,1,1,1],lambda);
    fluorophore.gaussian_405=fit_405([0,0,res_405(3),0,0,0,res_405(7),res_405(8),1,1,1,1,1,1],lambda);
    fluorophore.lipo_385=fit_385([0,0,0,res_385(4),0,0,1,1,res_385(9),res_385(10),1,1,1,1],lambda);
    fluorophore.lipo_405=fit_405([0,0,0,res_405(4),0,0,1,1,res_405(9),res_405(10),1,1,1,1],lambda);
    fluorophore.PpIX_636_385=fit_385([0,0,0,0,res_385(5),0,1,1,1,1,res_385(11),res_385(12),1,1],lambda);
    fluorophore.PpIX_636_405=fit_405([0,0,0,0,res_405(5),0,1,1,1,1,res_405(11),res_405(12),1,1],lambda);
    fluorophore.PpIX_620_385=fit_385([0,0,0,0,0,res_385(6),1,1,1,1,1,1,res_385(13),res_385(14)],lambda);
    fluorophore.PpIX_620_405=fit_405([0,0,0,0,0,res_405(6),1,1,1,1,1,1,res_405(13),res_405(14)],lambda);
    fluorophore.StO2=x_short(5);
    fluorophore.c_blood=x_short(4);
end