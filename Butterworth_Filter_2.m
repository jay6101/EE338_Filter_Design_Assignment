%Butterworth Analog LPF parameters
Wc = 1.074;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + 1i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - 1i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + 1i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - 1i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + 1i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - 1i*Wc*sin(pi/2 + pi/16+3*pi/8);

%Band Edge speifications
F_pl = 38.7e3;
F_sl = 42.7e3;
F_sh = 62.7e3;
F_ph = 66.7e3;

%Transformed Band Edge specs using Bilinear Transformation
F_samp = 260e3;         
W_pl = tan(F_pl/F_samp*pi);
W_sl = tan(F_sl/F_samp*pi); 
W_sh = tan(F_sh/F_samp*pi);
W_ph = tan(F_ph/F_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(W_pl*W_ph);
B = W_ph-W_pl;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, F_samp);
plot(f,abs(H))

yline(0.15,'-o','Magnitude = 0.15');
yline(0.85,'-o','Magnitude = 0.85');
yline(1,'-o','Magnitude = 1');
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