%% Processing of fluorescence data
clear all
close all
addpath('D:\Lyon thèse\soft\soft manip Arthur\functions')
clc
% cd 'C:\Users\auzel.NB411\Documents\Lyon thèse\soft\soft manip Arthur'
% dataDir='D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data';
% dataDir='D:\Lyon thèse\data\Manip fantomes\Phantom_20250205_NADH_FAD_TiO2_AR1_AR14\System_008_fluo';
% % dataDir='D:\Lyon thèse\data\Manip fantomes\2025_02_13_spectralon_mirror\System_008_fluo';
% dataDir='D:\Lyon thèse\data\Manip cochons Lyon\2025_02_21\System_008_fluo';
% dataDir='D:\Lyon thèse\data\Manip fantomes\2025_02_28_IGL2_pergen\System_008_fluo';
% dataDir='D:\Lyon thèse\soft\soft manip Arthur\correction_optical_p roperties\experimental_data';
% dataDir='D:\Lyon thèse\data\Manip fantomes\Phantom_20250205_NADH_FAD_TiO2_AR1_AR14\System_008_fluo\data_to_treat';
% dataDir="D:\Lyon thèse\data\Manip fantomes\2025_08_10_Bound_free_NADH";
% dataDir='D:\Lyon thèse\data\Manip fantomes\test angele yeast';
% dataDir='D:\Lyon thèse\data\Manip fantomes\Phantom_20241108_old_NADH\System_008_fluo';
%     dataDir='D:\Lyon thèse\data\Manip souris Marseille\2024_04_29\System_008_fluo\modified_data';
%     dataDir='D:\Lyon thèse\data\Manip fantomes\2024_07_05_AR1_AR14_FAD_TIO2_IL\System_008_fluo';
%     dataDir='D:\Lyon thèse\data\Manip fantomes\Phantom_20250205_NADH_FAD_TiO2_AR1_AR14\System_008_fluo\';
%     dataDir='D:\Lyon thèse\data\HYROX\2025_04_07_cells_PBS_fluo\System_008_fluo';
% dataDir="D:\Lyon thèse\data\Manip souris Marseille\2023_12_12\System_008_fluo";
%     dataDir='D:\Lyon thèse\soft\soft manip Arthur\correction_optical_properties\experimental_data';
    dataDir='D:\Lyon thèse\data\Manip cochons Lyon\2024_07_11\System_008_fluo';
% dataDir= "D:\Lyon thèse\data\Manip fantomes\2024_09_12_retinoic_acid\System_008_fluo";

    [data,path] = uigetfile('*fluo.mat',...
   'Select One or More Files', ...
   'MultiSelect', 'on',(dataDir));
data=convertCharsToStrings(data);

S385total=[];
S405total=[];
lambda=[];
raw_lambda=[];


figure
for k=1:length(data)
    load(fullfile(path,data(k)));
    LEDallummee=donnees_acq_brut.data(:,4:6);
    Power=donnees_acq_brut.data(:,7:7);
    idx385=find(LEDallummee(:,1)==1);
    idx405=find(LEDallummee(:,2)==1);
    Power385=Power(idx385(1));
    Power405=Power(idx405(1));
    S385=(signal_brut(idx385,:)-signal_brut(idx385-1,:))/Power385;
    S405=(signal_brut(idx405,:)-signal_brut(idx405-1,:))/Power405;
%     S385=signal_brut(idx385,:);
%     S405=signal_brut(idx405,:);
    S385total=[S385total;sum(S385(:,100:800))];
    S405total=[S405total;sum(S405(:,100:800))];
    if ~isempty(raw_lambda)
        lambda=raw_lambda;
    end
    subplot(211)
    hold on
    plot(lambda(100:800),S385(:,100:800))
    subplot(212)
    hold on
    plot(lambda(100:800),S405(:,100:800))
end
legend()


colors=jet(length(S385total(:,1)))*0.9;

figure('Name','Fluorescence laser')
legende=strrep(data,"_"," ");

figure

plot(lambda(100:800),S385total,'LineWidth',2)
legend(erase(legende,".mat"))
xlabel("Wavelength (nm)")
ylabel("Fluorescence intensity (a.u)")
set(gca,'Fontsize',20,'Fontweight','bold')


figure
for k=1:length(S385total(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda(100:800),S385total(k,:)/max(S385total(k,:)))
    legend(erase(legende,".mat"))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
    subplot(2,1,2)
    hold on
    plot(lambda(100:800),S405total(k,:)/max(S405total(k,:)))
%     legende=eraseBetween(legende," ",".mat");
    legend(erase(legende,["treated",".mat"]))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
end

figure
for k=1:length(S385total(:,1)) 
    plot(lambda(100:800),S405total(k,:))
    hold on
    plot(lambda(100:800),S405total(k,:)/3)
%     legende=eraseBetween(legende," ",".mat");
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
end
legend()
figure

for k=1:length(S385total(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda(100:800),S385total(k,:))
    legend(erase(legende,".mat"))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
    subplot(2,1,2)
    hold on
    plot(lambda(100:800),S405total(k,:))
%     legende=eraseBetween(legende," ",".mat");
    legend(erase(legende,["treated",".mat"]))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
end

%%
figure('Name','Fluorescence laser','Position',[100 100 1100 700])


for k=1:length(S385total(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda(100:800),S385total(k,:)./sum(S385total(k,:)),'LineWidth',2)
    legend(erase(legende,".mat"))
    xlabel("Wavelength (nm)")
    ylabel("Normalized F.I (a.u)")
    xlim([450 700])
    title('\it Excitation at 375 nm')
    set(gca,'FontSize',14,'FontWeight','bold')
    subplot(2,1,2)
    hold on
    plot(lambda(100:800),S405total(k,:)./sum(S405total(k,:)),'LineWidth',2)
%     legende=eraseBetween(legende," ",".mat");
    legend(erase(legende,["treated",".mat"]))
    xlabel("Wavelength (nm)")
    ylabel("Normalized F.I (a.u)")
    xlim([450 700])
    title('\it Excitation at 405 nm')
    set(gca,'FontSize',14,'FontWeight','bold')

end


figure('Name','Fluorescence laser')


for k=1:length(S385total(:,1)) 
    subplot(2,1,1)
    title('laser 385')
    hold on
    plot(lambda(100:800),S385total(k,:)./sum(S385total(k,:)))
    legend(erase(legende,".mat"))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
    set(gca,'FontSize',20)
    subplot(2,1,2)
    title('laser 405')
    hold on
    plot(lambda(100:800),S405total(k,:)./sum(S405total(k,:)))
%     legende=eraseBetween(legende," ",".mat");
    legend(erase(legende,["treated",".mat"]))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
    set(gca,'FontSize',20)
end


figure
plot(max(S385total,[],2)./max(S405total,[],2),'*')
legend(data)


%%
figure
for k=1:length(S385total(:,1)) 
    title('laser 385')
    hold on
    plot(lambda(100:800),S385total(k,:)./sum(S385total(k,:)),'LineWidth',2)
    legend(erase(legende,".mat"))
    xlabel("Wavelength (nm)")
    ylabel("Fluorescence intensity (a.u)")
    set(gca,'FontSize',30,'FontWeight','bold')
end



%% Processing of white light data
clear all 
close all
addpath('D:\Lyon thèse\soft\soft manip Arthur\functions')

file_spectralon_theo='D:\Lyon thèse\data\theoretical data\spectralon\Reflectance_values_array.txt';

% dataDir='D:\Lyon thèse\data\HYROX';
    dataDir='D:\Lyon thèse\data\Manip cochons Lyon\2023_04_04\System_008_fluo';
min_wl=480;
max_wl=645;

% 

% % dataDir="D:\Lyon thèse\data\Manip souris Marseille\2024_04_29\System_008_fluo";
legende=[];
name_spectralon="spectralon_000_diffuseReflectances.mat";
% we process data acquired after obtaining the source spectrum
% files_WL=name_WL(dataDir);
[data,path] = uigetfile('*diffuseReflectances.mat',...
   'Select One or More Files', ...
   'MultiSelect', 'on',(dataDir));

data=convertCharsToStrings(data);
WLShortTotal=[];
WLLongTotal=[];



for k=1:length(data)
%     read=contains(files_WL,data(k));
%     files_WL(read);
    WLight=load(fullfile(path,data(k)));
    WLallummee=WLight.donnees_acq_brut.data(:,4:6);
    idxLong=find(WLallummee(:,2)==1);
    idxShort=find(WLallummee(:,1)==1);
    WLShort=WLight.signal_brut(idxShort,:)-WLight.signal_brut(idxShort-1,:);
    WLLong=WLight.signal_brut(idxLong,:)-WLight.signal_brut(idxLong-1,:);
    WLLongTotal=[WLLongTotal;mean(WLLong)];
    WLShortTotal=[WLShortTotal;mean(WLShort)]; % we remove the last spectrum since it is the BG here, and we remove the last value of each spectrum which is not data

end
% Division by the source spectrum
if isfield(WLight,'lambda') 
    lambda=WLight.lambda;
else 
    lambda=WLight.raw_lambda;
end
% lambda=lambda(100:800);
legende=strrep(data,"_"," ");

figure
plot(lambda(100:800),WLShort(:,100:800))

%%
figure('Position', [100, 100, 800, 600]);

plot(lambda,WLShortTotal/max(WLShortTotal),'LineWidth',4)
legend(erase(legende,"diffuseReflectances.mat"))
grid on
xlabel("Wavelength (nm)")
ylabel("Reflectance")
set(gca,'FontSize',16,'FontWeight','bold')
xlim([435 750])

figure('Name','White Light DRS');
subplot(2,1,1)
hold on
plot(lambda,WLShortTotal(:,100:800))
legend(erase(legende,"diffuseReflectances.mat"))
xlabel("Wavelength (nm)")
ylabel("Reflectance")
subplot(2,1,2)
hold on
plot(lambda,WLLongTotal(:,100:800))
legend(erase(legende,"diffuseReflectances.mat"))
xlabel("Wavelength (nm)")
ylabel("Reflectance")

figure('Name','White Light DRS');
subplot(2,1,1)
hold on
plot(lambda,WLShortTotal(:,100:800)./sum(WLShortTotal(:,100:800),2))
legend(erase(legende,"diffuseReflectances.mat"))
xlabel("Wavelength (nm)")
ylabel("Reflectance")
subplot(2,1,2)
hold on
plot(lambda,WLLongTotal(:,100:800)./sum(WLLongTotal(:,100:800),2))
legend(erase(legende,"diffuseReflectances.mat"))
xlabel("Wavelength (nm)")
ylabel("Reflectance")


WLShortTotalReal=[];
WLLongTotalReal=[];
WLShortTotalReal2=[];
WLLongTotalReal2=[];
for i=1:length(WLShortTotal(:,1))
    WLnewShort=WLShortTotal(i,100:800);
    WLnewLong=WLLongTotal(i,100:800);
    if isfield(WLight,'name_spectralon')
        file_spectralon_exp=fullfile(dataDir,WLight.name_spectralon);
    else 
        file_spectralon_exp=fullfile(dataDir,name_spectralon);
    end
    [SsourceShort,SsourceLong]=calibration_spectralon(file_spectralon_exp,file_spectralon_theo);
    WLnewShort= WLnewShort./SsourceShort(100:800);
    WLnewLong=WLnewLong./SsourceLong(100:800);
    WLnewShort2= WLnewShort/max(WLnewShort);
    WLnewLong2=WLnewLong/max(WLnewLong);
    WLShortTotalReal2=[WLShortTotalReal2 ; WLnewShort2];
    WLLongTotalReal2=[WLLongTotalReal2 ; WLnewLong2];
    WLShortTotalReal=[WLShortTotalReal ; WLnewShort];
    WLLongTotalReal=[WLLongTotalReal ; WLnewLong];

end



figure('Name','White Light DRS');
plot(lambda,WLShortTotalReal,'LineWidth',2)
set(gca,'FontSize',20,'FontWeight','bold')
legend()
xlabel('Wavelength (nm)')
ylabel('Experimental reflectance (a.u)')
colors=jet(length(WLShortTotal(:,1)))*0.9;

figure


for k=1:length(WLShortTotal(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda,WLShortTotalReal(k,:),'color',colors(k,:))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance")
    subplot(2,1,2)
    hold on
    plot(lambda,WLLongTotalReal(k,:),'color',colors(k,:))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance")
end

figure('Name','White Light DRS with different perfusion time with oxygenated blood long distance');

for k=1:length(WLShortTotalReal2(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda,WLShortTotalReal2(k,:))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance normalized")
    subplot(2,1,2)
    hold on
    plot(lambda,WLLongTotalReal2(k,:))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance normalized")
end

figure
for k=1:length(WLShortTotalReal2(:,1)) 
    subplot(2,1,1)
    hold on
    plot(lambda,WLShortTotalReal(k,:)/sum(WLShortTotalReal(k,:)))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance normalized")
    subplot(2,1,2)
    hold on
    plot(lambda,WLLongTotalReal2(k,:)/sum(WLLongTotalReal(k,:)))
    legend(erase(legende,"diffuseReflectances.mat"))
    xlabel("Wavelength (nm)")
    ylabel("Reflectance normalized")
end