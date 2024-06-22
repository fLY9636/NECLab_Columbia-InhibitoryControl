function Filtered = Filter(Raw,Fs,order,fc,type)
[b, a] = butter(order,fc*2/Fs,type);
Filtered = filtfilt(b,a,Raw);
end
