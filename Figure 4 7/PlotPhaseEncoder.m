e1.MetaData = MetaData_All;
i = 56;
S2P_IDX = e1.MetaData.Sensor_time{i};
T = (1:length(S2P_IDX))/120;
i_NE = Filter(e1.MetaData.NE_470{i}, 120, 2, [0.4 0.8], 'bandpass');
% i_NE = Filter(e1.MetaData.NE_470{250}, 120, 2, [0.5 1.2], 'bandpass'); i_NE = zscore(i_NE(S2P_IDX));
i_Ach = Filter(e1.MetaData.Ach_470{i}, 120, 2, [0.4 0.8], 'bandpass');
% i_Ach = Filter(e1.MetaData.Ach_470{250}, 120, 2, [0.5 1.2], 'bandpass'); i_Ach = zscore(i_Ach(S2P_IDX)); 
figure; hold on
st = 1405.5; ed = 1420.5;
subplot(5, 1, 1); hold on
plot(T, i_NE, 'b'); plot(T, i_Ach, 'r'); xlim([st ed])
title('[0.4 0.8] BPF')
subplot(5, 1, 2); hold on
plot(T, angle(hilbert(i_NE)), 'b'); plot(T, angle(hilbert(i_Ach)), 'r'); xlim([st ed])
title('Hilbert phase angle BPF')
subplot(5, 1, 3); hold on
plot(T, 1-sin(abs(angle(hilbert(i_NE))-angle(hilbert(i_Ach)))/2), 'k'); xlim([st ed])
title('Phase synchrony')
subplot(5, 1, 4); hold on
plot(T, angle(hilbert(1-sin(abs(angle(hilbert(i_NE))-angle(hilbert(i_Ach)))/2))), 'k'); xlim([st ed])
title('Hilbert phase angle synchrony')
subplot(5, 1, 5); hold on
sync = 1-sin(abs(angle(hilbert(i_NE))-angle(hilbert(i_Ach)))/2);
sync = Filter(sync, 120, 2, 4, 'low');
hilbert_sync = [diff(angle(hilbert(zscore(sync)))); 0];
sync_hilbert = angle(hilbert(zscore(sync)));
spikes = NaN(size(sync)); 
spikes(find(hilbert_sync<=-5)) = 1; 
stem(T, spikes, '-k', 'Marker', 'none'); ylim([-0.5 2]); xlim([st ed])
title('Spikes')