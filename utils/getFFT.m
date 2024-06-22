function [P, PW, pxx, f, f1] = getFFT(input, fs, N, fig)
input = input(1:N);
P = abs(fft(input, length(input)));
[pxx, f1] = pspectrum(input, fs, 'power');
PW = pow2db(pxx);
f = (fs/length(input)):(fs/length(input)):(fs/2);
if strcmp(fig, 'on')==1
    figure
    subplot(2, 2, 1)
    grid on
    plot(f, P(1:(length(input)/2)))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title('FFT')
    subplot(2, 2, 3)
    xlabel('Logarithm frequency (Hz)')
    ylabel('Amplitude')
    grid on
    semilogx(f, P(1:(length(input)/2)))
    subplot(2, 2, [2 4])
    grid on
    plot(f1, PW)
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    title('Power Spectrum')
else
    ;
end
return