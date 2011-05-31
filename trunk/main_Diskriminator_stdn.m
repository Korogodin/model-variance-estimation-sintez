%clear all
clc
close all

Tc = 0.001; % Интервал накопления в корреляторе
Tk = 0.1; % Период такта работы дискриминаторов
K = fix(Tk / Tc); 
stdn_IQ = 600;

Nstat = 10000; 
qcno_dB = 33;
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ, Tc); % Амплитуда, соответсвующая выбранному qcno_dB

phase = rand(1,K)*2*pi;

d_stdn2_IQ = 0.6*stdn_IQ.^2/50;
SQ_stdn_IQ_extr = stdn_IQ.^2 + d_stdn2_IQ*(-50:100); % Точки экстраполяции stdn_IQ для дискриминатора
size_extr = length(SQ_stdn_IQ_extr);

Ud = zeros(1, size_extr);
sumI10 = nan(1, size_extr);
for nstat = 1:Nstat
    n_I = randn(1,K);
    n_Q = randn(1,K);
%     I = A_IQ .* cos(phase) + n_I.*stdn_IQ;
%     Q = - A_IQ .* sin(phase) + n_Q.*stdn_IQ;
    X2_arr = (A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2; % (I^2 + Q^2)_k
%     sumX = sum( sqrt(X2_arr) );
    sumX2 = sum( X2_arr );
%     sumX4 = sum( ((A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2).^2 );
    SQ_A_IQ_extr = A_IQ.^2 * (1 + 0.00*randn(1,1)); % Используемая экстраполяция A^2
    for nextr = 1:size_extr
        sumI10(nextr) = sum( besseli(1, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr(nextr)) ...
                    ./besseli(0, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr(nextr))  ...
                    .* sqrt(X2_arr) );
    end
    Ud_k = - K./SQ_stdn_IQ_extr + 0.5*(sumX2 + K*SQ_A_IQ_extr)./SQ_stdn_IQ_extr.^2  - sumI10 * sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr.^2;
    Ud = Ud + Ud_k;
    if ~mod(100*nstat, Nstat)
        fprintf('Complete %.0f percents \n', (100*nstat/Nstat));
    end
end
Ud = Ud / nstat;

hF = 0;
hF = figure(hF+1);
plot( 10*log10(stdn_IQ.^2) - 10*log10(SQ_stdn_IQ_extr), Ud);
xlabel('\delta\sigma^2_{IQ}, dB');
ylabel('U_{d,\sigma}');
grid on

