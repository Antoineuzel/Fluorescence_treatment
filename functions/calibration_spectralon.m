function [SsourceShort,SsourceLong]=calibration_spectralon(file_spectralon_exp,file_spectralon_theo)
SpectralonExp=load(file_spectralon_exp);
WLallummee=SpectralonExp.donnees_acq_brut.data(:,4:6);
idxLong=find(WLallummee(:,2)==1);
idxShort=find(WLallummee(:,1)==1);
SpecexpShort=SpectralonExp.signal_brut(idxShort,:)-SpectralonExp.signal_brut(idxShort-1,:);
SpecexpLong=SpectralonExp.signal_brut(idxLong,:)-SpectralonExp.signal_brut(idxLong-1,:);
%valeurs théoriques du spectralon
spectralon=importdata(file_spectralon_theo);
spectralon=spectralon.data;
if ~isfield(SpectralonExp,'lambda') 
    SpectralonExp.lambda=SpectralonExp.raw_lambda;
end
spectralonTheorique=interp1(spectralon(:,1),spectralon(:,2),SpectralonExp.lambda);
% Spectre_exp_norm=SpecexpShort(1,:)/max(SpecexpShort(1,:));
SsourceShort=sum(SpecexpShort)./(spectralonTheorique*length(idxShort));
SsourceLong=sum(SpecexpLong)./(spectralonTheorique*length(idxLong));
% Ssource=Ssource./max(Ssource);
% figure
% subplot(2,1,1)
% plot(SpectralonExp.lambda,SpecexpShort,SpectralonExp.lambda,SpecexpLong)
% ylabel("Source spectrum from both fibers")
% xlabel("Wavelength (nm)")
% subplot(2,1,2)
% plot(SpectralonExp.lambda,Ssource)
% ylabel("Normalized source spectrum")
% xlabel("Wavelength (nm)")
end