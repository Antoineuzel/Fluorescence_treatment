function [WLShort,WLLong,lambda,WLight,idx_min,idx_max,WLshort_all,WLlong_all] = load_reflectance(path,wl_min,wl_max)
%LOAD_REFLECTANCE Loads diffuse reflectance data from a *diffuseReflectances.mat file.
%   [WLShort,WLLong,lambda,WLight,idx_min,idx_max] = load_reflectance(path,wl_min,wl_max)
%   returns the mean raw signal (across acquisitions) for the short and
%   long wavelength LEDs, background-subtracted.
%   [...,WLshort_all,WLlong_all] additionally returns every individual
%   acquisition (rows), useful to inspect acquisition-to-acquisition
%   variability (see explore_single_measurement.m).
WLight=load(path);
if isfield(WLight,'lambda')
    lambda=WLight.lambda;
else
    lambda=WLight.raw_lambda;
end
idx_min=find(min((lambda-wl_min).^2)==(lambda-wl_min).^2);
idx_max=find(min((lambda-wl_max).^2)==(lambda-wl_max).^2);
lambda=lambda(idx_min:idx_max);
WLallummee=WLight.donnees_acq_brut.data(:,4:6);
idxLong=find(WLallummee(:,2)==1);
idxShort=find(WLallummee(:,1)==1);
WLshort_all=WLight.signal_brut(idxShort,idx_min:idx_max)-WLight.signal_brut(idxShort-1,idx_min:idx_max);
WLlong_all=WLight.signal_brut(idxLong,idx_min:idx_max)-WLight.signal_brut(idxLong-1,idx_min:idx_max);
WLShort=mean(WLshort_all,1);
WLLong=mean(WLlong_all,1);
end