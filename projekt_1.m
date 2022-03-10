clc;
close all;
clear all;

%% bod 1

[n,Wp] = cheb1ord(1000*2*pi,250*2*pi, 3, 20,'s');
%n - minimalni rad systemu
%Wp - zlomova frekvence

%wp - frekvence na obrazku = 1 kHz/rad*s^-1
%ws - frekvence na obrazku = 250 kHz/rad*s^-1
%|FP| = amplituda na obrazku = -3 dB
%|FS| = amplituda na obrazku = -20dB
%|F(nekonecno)| = 0dB

% je k tomu obrazek

% prototyb Čebyšev - proto funkce cheb1ord

%wp v normalnich jednotkach

[z,p,k] = cheb1ap(n,3)% prenos, trojka je tam Fp

k = k*10^(3/20); % zmenime gain, abychom posunuli primku nahoru - o 3 decibely

[b,a] = zp2tf(z,p,k)

[bt,at] = lp2hp(b,a,Wp * 1.1685) % posuneme do strany

[b,a] = zp2tf(z,p,k);

F1 = tf(bt,at) % vypoceteny prenos

% bode(F)
%grid on;


%% bod 2

C = 33.2 * 10^-9 ;

[num, den] = tfdata(F1);

k0 = num{1}(1);

k1 = den{1}(3);

k2 = den{1}(2);

R1_v = k2 / (k1* 2 * C);    % vypoctene hodnoty

R2_v = 2/(k2 * C);

R3 = 0;

R4 = inf;




%% bod 3

C = 33.2 * 10^-9 ;
sym R0;
sym w0;

w0 = 1000 * 2 * pi; %mezni frekvence pro RC clanek






C0 = 33.63 * 10^-9; % skutecne hdontoy
C1 = 33.97 * 10^-9;
C2 = 33.45 * 10^-9;

R0 = 4.7251 * 10^3;  %skutecne hdonoty
R1 = 1.284 * 10^3;
R2 = 8.799 * 10^3;

R0_v = 1/(w0*C0); % vypoctena hodnota

odchylka_R1 = (R1_v-R1)/R1; %overeni, jestli jsou namerene hodnoty malo odlisne od skutecnych
odchylka_R2 = (R2_v-R2)/R2;

num2 = [C1*C2*1, 0 , 0]
den2 = [C1*C2, (C1+C2)/R2 - 0 , 1/(R1*R2)]

F2 = tf(num2,den2) % prenos z namerenych hodnot

figure;
bode(F1) %vypoctene
hold on
bode(F2) %namerene
title('Bodeho charakteristika bez RC Clanku')

%% Prenos s RC clankem

odchylka_R0 = (R0_v-R0)/R0;


numRC_v = [1/(R0_v*C0)  0];
denumRC_v = [1/(R0_v*C0) 1]
F1_RC = tf(numRC_v, denumRC_v); % z vypoctenych hodnot
F1_celkovy = F1_RC*F1
bode(F1_celkovy) % celkovy prenos

numRC = [1/(R0*C0)  0];
denumRC = [1/(R0*C0) 1];

F2_RC = tf(numRC, denumRC); % z namerenych nodnot
F2_celkovy = F2_RC*F2;

figure;
bode(F2_celkovy) % celkovy prenos
hold on
title('Bodeho charakterstika s RC clankem')


f = [50 100 200 500 800 900 1000 1100 1200 1500 2000 5000 10000 20000]
f = f*2*pi;
U2 = [32*10^-3 36*10^-3 56*10^-3 272*10^-3 880*10^-3 1.1 1.35 1.67 1.99 2.73 2.73 2.13 2.05 2.05]
U2_celkovy = [0*10^-3 0*10^-3 23*10^-3 93*10^-3 280*10^-3 370*10^-3 480*10^-3 680*10^-3 800*10^-3 1.16 1.89 2.13 2.05 2.01]

F = U2/2; %2V
F_celkovy = U2_celkovy/2; %2V

% chceme vykreslit v logaritmitckych souradnicich


F = log10(F)*20;
F_celkovy = log10(F_celkovy)*20;

figure;
semilogx(f,F)
hold on
bode(F2) % prenos bez RC clanku se skutecnymi hodnotami
title('Bodeho charakterstika bez RC clanku')



figure;
semilogx(f,F_celkovy)
hold on
bode(F2_celkovy) % celkovy prenos se skutecnymi hodnotami
title('Bodeho charakterstika s RC clankem')

%% Prace s xlsx soubory



