%clear all
for p = 1:1
clc
close all

Tc = 0.001; % Интервал накопления в корреляторе
Tk = 1; % Период такта работы дискриминаторов
K = fix(Tk / Tc); 
stdn_IQ = 1.5;

Nstat = 20000; 
qcno_dB = 5+p;
[A_IQ, qcno] = qcno_change(qcno_dB, stdn_IQ, Tc); % Амплитуда, соответсвующая выбранному qcno_dB

qcno_dB_extr = qcno_dB; 
[A_IQ_extr, qcno] = qcno_change(qcno_dB_extr, stdn_IQ, Tc); % Точки экстраполяции А для дискриминатора

xA = nan(1,Nstat+1);
SQ_A_IQ_extr = A_IQ_extr.^2; 
xA(1) = 0+10*log10(SQ_A_IQ_extr);
SQ_A_IQ_extr = 10^(xA(1)/10);

xS = nan(1,Nstat+1);
SQ_stdn_IQ_extr = (stdn_IQ).^2 * (1 + 0.00*randn(1,1)); % Используемая экстраполяция СКО^2
xS(1) = 0+SQ_stdn_IQ_extr;
SQ_stdn_IQ_extr = xS(1);
for nstat = 1:Nstat
    phase = rand(1,K)*2*pi;
    n_I = randn(1,K);
    n_Q = randn(1,K);
%     I = A_IQ .* cos(phase) + n_I.*stdn_IQ;
%     Q = - A_IQ .* sin(phase) + n_Q.*stdn_IQ;
    X2_arr = (A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2; % (I^2 + Q^2)_k
    sumX2 = sum( X2_arr );
%     sumX = sum( sqrt(X2_arr) );
%     sumX2 = sum( X2_arr );
%     sumX4 = sum( ((A_IQ .* cos(phase) + n_I.*stdn_IQ).^2 + (- A_IQ .* sin(phase) + n_Q.*stdn_IQ).^2).^2 );
    
        sumI10 = sum( besseli(1, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr) ...
                    ./besseli(0, sqrt(X2_arr).*sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr)  ...
                    .* sqrt(X2_arr) );
    Ud_A = - K  + sumI10 ./ sqrt(SQ_A_IQ_extr);
    Ka = 0.004;
    xA(nstat+1) = xA(nstat) + Ka*Ud_A;
    SQ_A_IQ_extr = 10^(xA(nstat+1)/10);
    
    Ud_S = - K./SQ_stdn_IQ_extr + 0.5*(sumX2 + K*SQ_A_IQ_extr)./SQ_stdn_IQ_extr.^2  - sumI10 * sqrt(SQ_A_IQ_extr)./SQ_stdn_IQ_extr.^2;    
    if (qcno_dB > 20)
        Ks = 0.001;
    else
        Ks = 0.00001;
    end
    xS(nstat+1) = xS(nstat) + Ks*Ud_S;
    SQ_stdn_IQ_extr = xS(nstat+1);
    
    if ~mod(100*nstat, Nstat)
        fprintf('Complete %.0f percents \n', (100*nstat/Nstat));
    end
end
%Ud = Ud / nstat;
scr_ErrSys;
end
% hF = 0;
% hF = figure(hF+1);
% plot(xA - 10*log10(2*xS*Tc));
% xlabel('t, s')
% ylabel('qcno, dBHz')
% grid on
% 
% 
% hF = figure(hF+1);
% plot(xA);
% xlabel('t, s')
% ylabel('xA, dB')
% grid on
% 
% hF = figure(hF+1);
% plot(xS);
% xlabel('t, s')
% ylabel('xS')
% grid on
