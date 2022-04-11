%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 4                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%Variables
dt = 1/1000; % Time step in s
time = linspace(0,1,1000);
Vout = zeros(1,length(time));

R1 = 1; R2 = 2; R4 = 0.1; R0 = 1000;
Cc = 0.25; L = 0.2; alpha = 100; beta = 100; gamma = 100; Vin = 1;
n = 5;
Cn = 10e-6; 
In = 0.001 * randn(1,length(time));

global G C b
C = sparse(n,n); G = sparse(n,n); b = sparse(n,1);

R3 = 20; % Temporary value for R3, my part 2 doesn't work yet

res(1,2,R1); res(2,0,R2); res(3,0,R3);  
res(4,5,R4); res(5,0,R0); cap(1,2,Cc);
cap(3,0,Cn); % Added noise capacitor
vol(1,0,Vin);
ind(2,3,L);
vcvs(4,0,3,0,alpha/R3+(beta/R3)^2+(gamma/R3)^3);

oldV = zeros(8,1); V = zeros(8,1);

for index = 1:length(time)
    cur(3,0,In(index)); % Adding noise current
    b(6) = exp(-(time(index)-0.06)^2/(2*0.03^2)); 
    
    oldV = V;
    A = C/dt+G;

    V = A\(b+C*oldV/dt);
    
    Vin(index) = V(1);   
    Vout(index) = V(5);
    
    plot(time, Vout,'g')
    title('Output Voltage With Noise Over Time')
    xlabel('Time (s)')
    ylabel('Output Voltage (V)')
    drawnow
end

FFTin = fft(Vin);
FFTout = fft(Vout);
FFTSin = fftshift(abs(FFTin));
FFTSout = fftshift(abs(FFTout));
x = linspace(-500,500,length(time));

figure
plot(x,FFTSin) 
hold on
plot(x,FFTSout) 
title('Input and Output Frequency Response')
legend('Input', 'Output')
xlabel('Frequency (f)')
ylabel('Magnitude (dB)')
