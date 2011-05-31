%p = 1;
q(p) = qcno_dB;

SKOhxA(p) = sqrt( sum( (xA - 20*log10(A_IQ)).^2 ) / (nstat+1) );
SKOhxS(p) = sqrt( sum( (xS - stdn_IQ.^2).^2 ) / (nstat+1) );

SKOqcnodB(p) = sqrt( sum( (xA - 10*log10(2*xS*Tc) - qcno_dB).^2 ) / (nstat+1) );


hF = 0;
hF = figure(hF+1);
plot(q, SKOhxA);
xlabel('q_{c/n0}, dBHz')
ylabel('MSE A_{IQ}^2 estimation, dB')
grid on


hF = figure(hF+1);
plot(q, SKOhxS);
xlabel('q_{c/n0}, dBHz')
ylabel('MSE \sigma_{IQ}^2 estimation')
grid on

hF = figure(hF+1);
plot(q, SKOqcnodB);
xlabel('q_{c/n0}, dBHz')
ylabel('MSE q_{c/n0,} estimation, dBHz')
grid on

drawnow;