clear all
clc
close all

settings;

Nmod = fix(Tmod / Tc); 

stdn_IQ_mean = (stdn_IQ_min + stdn_IQ_max)/2; % Среднее значение
stdn_IQ = stdn_IQ_min:d_stdn_IQ:stdn_IQ_max; % Массив значений stdn_IQ
Nstdn = size(stdn_IQ,2); % Число рассматриваемых stdn_IQ
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ_mean, Tc); % Амплитуда, соответсвующая выбранному qcno_dB

init_arrays;
phase = rand(1,Nmod)*2*pi;
A_IQ_izm = A_IQ*1.0;

for nmod = 1:Nmod
    n_I = randn(1,Nstdn);
    n_Q = randn(1,Nstdn);
    I = A_IQ * cos(phase(nmod)) + n_I.*stdn_IQ;
    Q = - A_IQ * sin(phase(nmod)) + n_Q.*stdn_IQ;
    X2 = X2 + I.^2 + Q.^2;    X4 = X4 + (I.^2 + Q.^2).^2;
end

SQ_stdn_izm_1 = (X2 + Nmod*A_IQ_izm.^2 + sqrt(X2.^2 - 6*Nmod*A_IQ_izm.^2.*X2 + Nmod^2*A_IQ_izm.^4))/ (4*Nmod);
SQ_stdn_izm_2 = (X2/Nmod - A_IQ_izm.^2)/2;
SQ_stdn_izm_3 = (X2/Nmod - sqrt(2*(X2/Nmod).^2 - X4/Nmod))/ (2);
SQ_stdn_izm_4 = (X2/Nmod + A_IQ_izm.^2 - sqrt(2*(X2/Nmod).^2 - X4/Nmod + A_IQ_izm.^4 - 2*A_IQ_izm.^2.*SQ_stdn_izm_3) )/ (2);

hF = 0;
hF = figure(hF+1);
plot(stdn_IQ, stdn_IQ, stdn_IQ, sqrt(SQ_stdn_izm_1), stdn_IQ, sqrt(SQ_stdn_izm_2), ...
    stdn_IQ, sqrt(SQ_stdn_izm_3), stdn_IQ, sqrt(SQ_stdn_izm_4));
legend('ist', '1','2','3', '4')
xlabel('stdn_{IQ}');
ylabel('sqrt SQ\_stdn\_izm');

hF = figure(hF+1);
plot(stdn_IQ, stdn_IQ.^2, stdn_IQ, SQ_stdn_izm_1, stdn_IQ, SQ_stdn_izm_2, ...
    stdn_IQ, SQ_stdn_izm_3, stdn_IQ, SQ_stdn_izm_4);
legend('ist', '1','2','3', '4')
xlabel('stdn_{IQ}');
ylabel('SQ\_stdn\_izm');