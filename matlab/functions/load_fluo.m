function [fluo_385_total,fluo_405_total,lambda,fluo385_all,fluo405_all] = load_fluo(path,wl_min,wl_max)
%LOAD_FLUO Loads fluorescence data from a *fluo.mat file.
%   [f385,f405,lambda] = load_fluo(path,wl_min,wl_max) returns the mean
%   spectrum (across acquisitions) for each laser, background-subtracted
%   and normalized by power*acquisition time.
%   [...,fluo385_all,fluo405_all] additionally returns every individual
%   acquisition (rows), useful to inspect acquisition-to-acquisition
%   variability (see explore_single_measurement.m).
    fluo=load(path);
    if isfield(fluo,'lambda')
    lambda=fluo.lambda;
    else
    lambda=fluo.raw_lambda;
    end
    idx_min=find(min((lambda-wl_min).^2)==(lambda-wl_min).^2);
    idx_max=find(min((lambda-wl_max).^2)==(lambda-wl_max).^2);
    lambda=lambda(idx_min:idx_max);
    LEDallummee=fluo.donnees_acq_brut.data(:,4:6);
    idx385=find(LEDallummee(:,1)==1);
    Power=fluo.donnees_acq_brut.data(:,7:7);
    idx405=find(LEDallummee(:,2)==1);
    Power385=Power(idx385(1));
    Power405=Power(idx405(1));
    time_acq=fluo.donnees_acq_brut.data(1,8)/1000;
    fluo385=(fluo.signal_brut(idx385,:)-fluo.signal_brut(idx385-1,:))/(Power385*time_acq);
    fluo405=(fluo.signal_brut(idx405,:)-fluo.signal_brut(idx405-1,:))/(Power405*time_acq);
    fluo385_all=fluo385(:,idx_min:idx_max);
    fluo405_all=fluo405(:,idx_min:idx_max);
    fluo_385_total=mean(fluo385_all,1);
    fluo_405_total=mean(fluo405_all,1);
end