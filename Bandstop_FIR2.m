F_samp = 260e3;

%Band Edge speifications
F_pl = 38.7e3;
F_sl = 42.7e3;
F_sh = 62.7e3;
F_ph = 66.7e3;

w_pl = F_pl*2*pi/F_samp;
w_sh = F_sh*2*pi/F_samp;
w_sl = F_sl*2*pi/F_samp;
w_ph = F_ph*2*pi/F_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end


N_min = ceil((A-7.95) / (2.285*(4*2*pi/260)));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 38;

%Ideal bandstop impulse response of length "n"
w_c1 = (w_pl + w_sl)/2;
w_c2 = (w_ph + w_sh)/2;
bs_ideal =  ideal_lp(pi,n) -ideal_lp(w_c2,n) + ideal_lp(w_c1,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, F_samp);
figure;
plot(f,abs(H));
xlim([20000,130000]);
yline(0.15,'-o','Magnitude = 0.15');
yline(0.85,'-o','Magnitude = 0.85');
yline(1.15,'-o','Magnitude = 1.15');
xline(38.7e3,'-b','f = 38.7kHz');
xline(42.7e3,'-b','f = 42.7kHz');
xline(62.7e3,'-b','f = 62.7kHz');
xline(66.7e3,'-b','f = 66.7kHz');
xlim([20000,120000]);
ylim([0,1.5]);
xlabel('Frequency (in Hz)');
ylabel('Magnitude');
title('Magnitude Plot');
hold on
grid on
figure;
plot(bs_ideal)
xlabel('Samples');
ylabel('Amplitude');
title('Time Domain');
hold on
grid on 