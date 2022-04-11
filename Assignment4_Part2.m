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
Vinput = zeros(1,length(time));

R1 = 1; R2 = 2; R4 = 0.1; R0 = 1000;
Cc = 0.25; L = 0.2; alpha = 100; Vin = 1;
n = 5;
global G C F
% C = sparse(n,n); G = sparse(n,n); F = sparse(n,1);
C = zeros(n,n); G = zeros(n,n); F = zeros(n,1);
R3 = 32.9; % Value from part 2

res(1,2,R1); res(2,0,R2); res(3,0,R3);  
res(4,5,R4); res(5,0,R0); cap(1,2,Cc);
vol(1,0,Vin);
ind(2,3,L);
vcvs(4,0,3,0,alpha/R3);

oldV = zeros(8,1); V = zeros(8,1);

for index = 1:length(time)
%     F(6) = PWL(time(index));
%     F(6) = sin(2*pi*10*time(index));
    F(6) = exp(-(time(index)-0.06)^2/(2*0.03^2));

    oldV = V;
    A = C/dt+G;

    V = A\(F+C*oldV/dt);
       
    Vout(index) = V(5);
    Vinput(index) = V(1);
    
    plot(time, Vinput,'b')
    title('Input Voltage Over Time')
    xlabel('time (s)')
    ylabel('Input Voltage (V)')
    hold on
    plot(time, Vout,'g')
    title('Output Voltage Over Time')
    xlabel('Time (s)')
    ylabel('Output Voltage (V)')
    drawnow
    hold off
end

FFTin = fft(Vinput);
FFTout = fft(Vout);
FFTSin = fftshift(abs(FFTin));
FFTSout = fftshift(abs(FFTout));
x = linspace(-500,500,1000);

figure
plot(x,FFTSin) 
hold on
plot(x,FFTSout) 
title('Input and Output Frequency Response')
legend('Input', 'Output')
xlabel('Frequency (f)')
ylabel('Magnitude (dB)')





