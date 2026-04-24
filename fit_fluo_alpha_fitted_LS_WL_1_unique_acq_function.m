clear all
close all
addpath("theoretical_data")
addpath('D:\Lyon thèse\soft\soft manip Arthur\functions')
dataDir='D:\Lyon thèse\data\Manip cochons Lyon\2024_07_17\System_008_fluo';
min_wl_reflectance=480;
max_wl_reflectance=645;
min_wl_fluo=480;
max_wl_fluo=645;
[FAD385total,FAD405total,lambda]=load_fluo("D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data\flavine_fluo.mat",min_wl_fluo,max_wl_fluo);
[NADH385total,NADH405total,lambda_NADH]=load_fluo("D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data\NADH_fluo_4.mat",min_wl_fluo,max_wl_fluo);
% [NADH385total,NADH405total,lambda_NADH]=load_fluo("D:\Lyon thèse\data\Manip fantomes\Phantom_20250205_NADH_FAD_TiO2_AR1_AR14\System_008_fluo\phantom13_NADH_after_003_fluo.mat",min_wl_fluo,max_wl_fluo);
%     dataDir='D:\Lyon thèse\data\Manip souris Marseille\2024_04_29\System_008_fluo\modified_data';
% dataDir="D:\Lyon thèse\data\Manip fantomes\2024_09_12_retinoic_acid\System_008_fluo";
fluorophores_385=[FAD385total;NADH385total];
fluorophores_405=[FAD405total;NADH405total];
% % 
[data,path] = uigetfile("*fluo.mat",...
   'Select One or More Files', ...
   'MultiSelect', 'on',(dataDir));
% 
% selectedFiles = uipickfiles('FilterSpec', fullfile(dataDir, '*fluo.mat'), ...
%                            'Prompt', 'Select One or More Files');

data=convertCharsToStrings(data);

% 

name_spectralon='spectralon_000_diffuseReflectances';
isSmallSpectralon='True';
data=convertCharsToStrings(data);
for i =1:length(data)
[S385total(i,:),S405total(i,:),S385_corrected(i,:),S405_corrected(i,:),fluorophore,lambda,res_385_correc_Kim_full_exp_tot(i,:),res_405_correc_Kim_full_exp_tot(i,:),residu] = corrected_fluo_LS_WL_1(dataDir,data(i),name_spectralon,isSmallSpectralon,min_wl_reflectance,max_wl_reflectance,min_wl_fluo,max_wl_fluo,fluorophores_385,fluorophores_405);
FAD_405(i,:)=fluorophore.flavine_405_exp;
FAD_385(i,:)=fluorophore.flavine_385_exp;
NADH_405(i,:)=fluorophore.NADH_405_exp;
NADH_385(i,:)=fluorophore.NADH_385_exp;
end
FAD_405_tot=sum(FAD_405,2);
FAD_385_tot=sum(FAD_385,2);
NADH_405_tot=sum(NADH_405,2);
NADH_385_tot=sum(NADH_385,2);
redox_405=sum(FAD_405,2)./(sum(FAD_405,2)+sum(NADH_405,2));
redox_385=sum(FAD_385,2)./(sum(FAD_385,2)+sum(NADH_385,2));


figure('Position',[100 100 800 600])

subplot(2,1,1)
plot(lambda,S385_corrected/sum(S385_corrected),'*',lambda,S385total/sum(S385total))
% xlim([min(lambda_study) sum(lambda_study)])
grid on
yl=ylim;
legend('Experimental data corrected','Experimental data')
xlabel('Wavelength (nm)')
ylabel('A.U')
title('laser 385')
subplot(2,1,2)
plot(lambda,S405_corrected/sum(S405_corrected),'*',lambda,S405total/sum(S405total))
yl=ylim;
grid on
% xlim([min(lambda_study) sum(lambda_study)])
legend('Experimental data corrected','Experimental data')
title('laser 405')
xlabel('Wavelength (nm)')
ylabel('A.U')


figure('Position', [100, 100, 800, 600]);
plot(lambda,S405total/sum(S405total),'*','Color',[0.7059, 0.5490, 0.4000], 'MarkerSize', 5)
hold on
plot(lambda,S405_corrected/sum(S405_corrected),'*','Color',[0, 0.4470, 0.7410],'MarkerSize', 5)
yl=ylim;
grid on
% xlim([min(lambda_study) sum(lambda_study)])
legend('Experimental data no correction','Experimental data corrected')
xlabel('Wavelength (nm)')
ylabel('Normalized fluorescence intensity (a.u)')
set(gca,'Fontsize',15,'FontWeight','bold')
xlim([min_wl_fluo,max_wl_fluo])
%%
sum_fluo_405=fluorophore.flavine_405_exp+fluorophore.gaussian_405+fluorophore.lipo_405+fluorophore.NADH_405_exp+fluorophore.PpIX_620_405+fluorophore.PpIX_636_405;
sum_fluo_385=fluorophore.flavine_385_exp+fluorophore.gaussian_385+fluorophore.lipo_385+fluorophore.NADH_385_exp+fluorophore.PpIX_620_385+fluorophore.PpIX_636_385;
diff=sum(abs(sum_fluo_405-S405_corrected));

figure('Position', [100, 100, 800, 600]);

p=plot(lambda,S405_corrected/sum(S405_corrected),'*',lambda,sum_fluo_405/sum(S405_corrected),lambda,fluorophore.NADH_405_exp/sum(S405_corrected),lambda,fluorophore.flavine_405_exp/sum(S405_corrected),lambda,fluorophore.gaussian_405/sum(S405_corrected),lambda,fluorophore.lipo_405/sum(S405_corrected),lambda,fluorophore.PpIX_636_405/sum(S405_corrected),lambda,fluorophore.PpIX_620_405/sum(S405_corrected));
% p=plot(lambda,S405_corrected,'*',lambda,sum_fluo_405,lambda,fluorophore.NADH_405_exp,lambda,fluorophore.flavine_405_exp,lambda,fluorophore.gaussian_405,lambda,fluorophore.lipo_405,lambda,fluorophore.PpIX_636_405,lambda,fluorophore.PpIX_620_405);
legend('Experimental data corrected','Sum of fluorophores','NADH','FAD','Protein bound FMN','Lipopigments','PpIX 636','PpIX 620')
xlabel('Wavelength (nm)')
ylabel('Normalized fluorescence intensity (a.u)')
xlim([min_wl_fluo,max_wl_fluo])

yl=ylim;
grid on
% xlim([min(lambda_study) sum(lambda_study)])
p(2).LineWidth=4;
p(3).LineWidth=4;
p(4).LineWidth=4;
p(5).LineWidth=4;
p(6).LineWidth=4;
p(7).LineWidth=4;
p(8).LineWidth=4;
set(gca,'FontSize',16,'FontWeight','bold')





figure('Position', [100, 100, 1100, 700]);
subplot(211)
title('\it Excitation at 375 nm')
hold on
p=plot(lambda,S385_corrected/sum(S385_corrected),'*',lambda,sum_fluo_385/sum(S385_corrected),lambda,fluorophore.NADH_385_exp/sum(S385_corrected),lambda,fluorophore.flavine_385_exp/sum(S385_corrected),lambda,fluorophore.gaussian_385/sum(S385_corrected),lambda,fluorophore.lipo_385/sum(S385_corrected),lambda,fluorophore.PpIX_636_385/sum(S385_corrected),lambda,fluorophore.PpIX_620_385/sum(S385_corrected));

% p=plot(lambda,S405_corrected,'*',lambda,sum_fluo_405,lambda,fluorophore.NADH_405_exp,lambda,fluorophore.flavine_405_exp,lambda,fluorophore.gaussian_405,lambda,fluorophore.lipo_405,lambda,fluorophore.PpIX_636_405,lambda,fluorophore.PpIX_620_405);
legend('Experimental data corrected','Sum of fluorophores','NADH','FAD','Protein bound FMN','Lipopigments','PpIX 636','PpIX 620')
xlabel('Wavelength (nm)')
ylabel('Normalized F.I (a.u)')
xlim([min_wl_fluo,max_wl_fluo])

yl=ylim;
grid on
% xlim([min(lambda_study) sum(lambda_study)])
p(2).LineWidth=4;
p(3).LineWidth=4;
p(4).LineWidth=4;
p(5).LineWidth=4;
p(6).LineWidth=4;
p(7).LineWidth=4;
p(8).LineWidth=4;
set(gca,'FontSize',14,'FontWeight','bold')

subplot(212)
title('\it Excitation at 405 nm')
hold on
p=plot(lambda,S405_corrected/sum(S405_corrected),'*',lambda,sum_fluo_405/sum(S405_corrected),lambda,fluorophore.NADH_405_exp/sum(S405_corrected),lambda,fluorophore.flavine_405_exp/sum(S405_corrected),lambda,fluorophore.gaussian_405/sum(S405_corrected),lambda,fluorophore.lipo_405/sum(S405_corrected),lambda,fluorophore.PpIX_636_405/sum(S405_corrected),lambda,fluorophore.PpIX_620_405/sum(S405_corrected));
    
% p=plot(lambda,S405_corrected,'*',lambda,sum_fluo_405,lambda,fluorophore.NADH_405_exp,lambda,fluorophore.flavine_405_exp,lambda,fluorophore.gaussian_405,lambda,fluorophore.lipo_405,lambda,fluorophore.PpIX_636_405,lambda,fluorophore.PpIX_620_405);
legend('Experimental data corrected','Sum of fluorophores','NADH','FAD','Protein bound FMN','Lipopigments','PpIX 636','PpIX 620')
xlabel('Wavelength (nm)')
ylabel('Normalized F.I (a.u)')
xlim([min_wl_fluo,max_wl_fluo])

yl=ylim;
grid on
% xlim([min(lambda_study) sum(lambda_study)])
p(2).LineWidth=4;
p(3).LineWidth=4;
p(4).LineWidth=4;
p(5).LineWidth=4;
p(6).LineWidth=4;
p(7).LineWidth=4;
p(8).LineWidth=4;
set(gca,'FontSize',14,'FontWeight','bold')

%%
figure('Position', [100, 100, 800, 600]);

p = plot(lambda, NADH385total/max(NADH385total), ...
         lambda, fluorophore.flavine_385_exp/max(fluorophore.flavine_385_exp), ...
         lambda, fluorophore.gaussian_385/max(fluorophore.gaussian_385), ...
         lambda, fluorophore.lipo_385/max(fluorophore.lipo_385), ...
         lambda, fluorophore.PpIX_636_385/max(fluorophore.PpIX_636_385), ...
         lambda, fluorophore.PpIX_620_385/max(fluorophore.PpIX_620_385));

% Ensuite appliquer les couleurs :
p(1).Color = [0.9290 0.6940 0.1250];  % Jaune
p(2).Color = [0.4940 0.1840 0.5560];  % Violet
p(3).Color = [0.4660 0.6740 0.1880];  % Vert
p(4).Color = [0.3010 0.7450 0.9330];  % Bleu clair
p(5).Color = [0.6350 0.0780 0.1840];  % Rouge foncé
p(6).Color = [0 0.4470 0.7410];       % Bleu

xlim([480,645])
p(1).LineWidth=6;
p(2).LineWidth=6;
p(3).LineWidth=6;
p(4).LineWidth=6;
p(5).LineWidth=6;
p(6).LineWidth=6;
legend('NADH','FAD','Protein bound FMN','Lipopigments','PpIX 636','PpIX 620')
grid on
xlabel('Wavelength (nm)')
ylabel('Normalized fluorescence intensity (a.u)')
set(gca,'FontSize',20,'FontWeight','bold')

% cd('D:\Lyon thèse\rédaction\manuscrit_final\Part5\chap3_content')
