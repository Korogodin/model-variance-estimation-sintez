% clear all
% clc
% close all

Tc = 0.001; % Correlation period

Tmod = 0.3; % Время моделирования

qcno_dB = 15; % q c/no in dBHz

stdn_IQ_min = 1; % Минимальное значение СКО квадратур
stdn_IQ_max = 2; % Максимальное --//--
d_stdn_IQ = 0.01;  % Шаг по stdn_IQ


Nmod = fix(Tmod / Tc); 

stdn_IQ_mean = (stdn_IQ_min + stdn_IQ_max)/2; % Среднее значение
stdn_IQ = stdn_IQ_min:d_stdn_IQ:stdn_IQ_max; % Массив значений stdn_IQ
stdn_IQ = stdn_IQ_mean;
Nstdn = size(stdn_IQ,2); % Число рассматриваемых stdn_IQ
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ_mean, Tc); % Амплитуда, соответсвующая выбранному qcno_dB

I = nan(1, 1);
Q = nan(1, 1);
X2 = zeros(1,1);
X4 = zeros(1,1);
phase = rand(1,Nmod)*2*pi;
A_IQ_izm = (A_IQ/8):(A_IQ/20):(4*A_IQ);
NA = size(A_IQ_izm,2);

I0 = ones(1,Nstdn);
for nmod = 1:Nmod
    n_I = randn(1,1);
    n_Q = randn(1,1);
    I = A_IQ * cos(phase(nmod)) + n_I.*stdn_IQ_mean;
    Q = - A_IQ * sin(phase(nmod)) + n_Q.*stdn_IQ_mean;
    X2 = X2 + I.^2 + Q.^2;    X4 = X4 + (I.^2 + Q.^2).^2;
    I0 = I0.*besseli(0, A_IQ_izm.*sqrt(I.^2 + Q.^2)./stdn_IQ.^2);
  
end

SQ_A_izm_1 = (X2/Nmod - 2*stdn_IQ^2);

K = Nmod;
p_aposter =  exp(-(X2 + K*A_IQ_izm.^2)/2./stdn_IQ.^2) .* I0;
p_aposter_s =  sum(p_aposter)*(A_IQ/20);
p_aposter = p_aposter / p_aposter_s;
hF = 0;
hF = figure(hF+1);
plot(A_IQ_izm, p_aposter, ...
        [sqrt(SQ_A_izm_1) sqrt(SQ_A_izm_1)], [0 max(p_aposter)]);
xlabel('A_{IQ,izm}')
ylabel('p( A_{IQ,izm}, \sigma_{IQ,ist}^2 | Y_1^K )')