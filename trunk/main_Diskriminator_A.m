%clear all
clc
close all

Tc = 0.001; % Интервал накопления в корреляторе
Tk = 0.1; % Период такта работы дискриминаторов
K = fix(Tk / Tc); 
stdn_IQ = 1.5;

Nstat = 10000; 
qcno_dB = 15;
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ, Tc); % Амплитуда, соответсвующая выбранному qcno_dB

phase = rand(1,K)*2*pi;

qcno_dB_extr = 1:0.5:70; 
[A_IQ_extr, qcno] = qcno_change(qcno_dB_extr, stdn_IQ, Tc); % Точки экстраполяции А для дискриминатора
SQ_A_IQ_extr = A_IQ_extr.^2; 
size_extr = length(SQ_A_IQ_extr);

Ud = zeros(1, size_extr);
sumI10 = nan(1, size_extr);
for nstat = 1:Nstat
    n_I = randn(1,K);
    n_Q = randn(1,K);
%     I = A_IQ .* cos(phase) + n_I.*stdn_IQ;
%     Q = - A_IQ .* sin(phase) + n_Q.*stdn_IQ;
    X2_arr = (A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2; % (I^2 + Q^2)_k
%     sumX = sum( sqrt(X2_arr) );
%     sumX2 = sum( X2_arr );
%     sumX4 = sum( ((A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2).^2 );
    SQ_stdn_IQ_extr = (stdn_IQ).^2 * (1 + 0.00*randn(1,1)); % Используемая экстраполяция СКО^2
    for nextr = 1:size_extr
        sumI10(nextr) = sum( besseli(1, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr(nextr))./SQ_stdn_IQ_extr) ...
                    ./besseli(0, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr(nextr))./SQ_stdn_IQ_extr)  ...
                    .* sqrt(X2_arr) );
    end
    Ud_k = - K  + sumI10 ./ sqrt(SQ_A_IQ_extr);
    Ud = Ud + Ud_k;
    if ~mod(100*nstat, Nstat)
        fprintf('Complete %.0f percents \n', (100*nstat/Nstat));
    end
end
Ud = Ud / nstat;

hF = 0;
hF = figure(hF+1);
plot(20*log10(A_IQ)-20*log10(A_IQ_extr), Ud);
xlabel('\delta A^2_{IQ}');
ylabel('Ud');
grid on

hF = figure(hF+1);
plot(A_IQ-A_IQ_extr, Ud);
xlabel('\delta A_{IQ}');
ylabel('Ud');
grid on

hF = figure(hF+1);
plot(20*log10(A_IQ_extr), Ud, 20*log10([A_IQ A_IQ]), [min(Ud) max(Ud)]);
xlabel('A_{IQ,extr}');
ylabel('Ud');
grid on