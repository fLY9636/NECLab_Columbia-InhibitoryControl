function [phi_NE, phi_Ach, spikes, sync, spikes_01, sync_hilbert] = CalcSpikes(NE, Ach, fs)
% Input: NE and Ach band-pass z-scored vectors
% Output: spikes, spikes_01(with NaN values replaced by 0s), sync_hilbert
phi_NE = angle(hilbert(NE));
phi_Ach = angle(hilbert(Ach));
phasediff = angle(hilbert(Ach))-angle(hilbert(NE));
sync = 1-sin(abs(angle(hilbert(Ach))-angle(hilbert(NE)))/2);
sync = Filter(sync, fs, 2, 4, 'low');
hilbert_sync = [diff(angle(hilbert(zscore(sync)))); 0];
sync_hilbert = angle(hilbert(zscore(sync)));
spikes = NaN(size(sync)); 
spikes(find(hilbert_sync<=-5)) = 1; 
spikes_01 = spikes;
spikes_01(find(isnan(spikes_01))) = 0;
return