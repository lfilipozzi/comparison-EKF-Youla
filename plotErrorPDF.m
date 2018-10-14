function plotErrorPDF(signal_est,signal_true,nb_points)
%PLOTERRORPDF Shows the PDF of the estimation error

% Resample true signal
signal_true = resample(signal_true,signal_est.Time);

% Compute estimation error
error   = signal_est.Data - signal_true.Data;
x_range = linspace(min(error),max(error),nb_points);

% Compute PDF
[f,x_pdf]  = hist(error,x_range);
plot(x_pdf,f/((max(x_pdf)-min(x_pdf))*sum(f)/length(x_pdf)));
xlabel('Estimation error');
ylabel('Probability density function');
end

