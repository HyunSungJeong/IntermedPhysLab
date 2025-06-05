function freq_3dB = plot_PreampGain(Gain, num_linfit)

    path = [fileparts(mfilename('fullpath')),filesep,'data',filesep,'A1',filesep,'A1Gain',sprintf('%d',Gain),'.txt'];

    data = importdata(path);
    if isstruct(data)
        data = data.data;
    end
    
    freq = data(:,1);
    Gain_dB = 10*log10(data(:,2));
    Gain_dB_low = Gain_dB - 10*log10(data(:,2) - data(:,3));
    Gain_dB_high = 10*log10(data(:,2) - data(:,3)) - Gain_dB;

    freq_lin = log10(freq(end-num_linfit+1 : end));
    Gain_lin = Gain_dB(end-num_linfit+1 : end);

    coeff = polyfit(freq_lin, Gain_lin, 1);
    freq_fit = linspace(freq_lin(1), freq_lin(end)+0.1, 10);
    Gain_fit = polyval(coeff, freq_fit);
    freq_fit = power(10, freq_fit);

    freq_3dB = (10*log10(Gain)-3-coeff(2))/coeff(1);
    freq_3dB = power(10,freq_3dB);

    Xmin = freq(1)*0.5;
    Xmax = freq(end)*7;
    Ymin = min(Gain_dB)-1;
    Ymax = max(Gain_dB)+1;

    figure;
    hold on;
    ax.XLabel.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.XAxis.MinorTick = 'on';
    set(gca, 'XScale', 'log');
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    
    legends = zeros(4,1);
    legendNames = {['$\mathrm{gain: \ ',sprintf('%.1f', Gain) ,' }$'], '$\mathrm{best \ fit}$', '$\mathrm{3dB \ frequency}$', '$\mathrm{measured \ gain}$'};
    legends(1) = plot(linspace(Xmin, Xmax, 10), repmat(10*log10(Gain), [1,10]), '--', 'linewidth', 1, 'Color', 'black');
    legends(2) = plot(freq_fit, Gain_fit, '--', 'linewidth', 1, 'Color', 'red');
    legends(3) = plot(repmat(freq_3dB, [1,10]), linspace(Ymin, 10*log10(Gain)-3, 10), '--', 'linewidth', 1, 'Color', 'blue');
    err = errorbar(freq, Gain_dB, Gain_dB_low, Gain_dB_high, 'o', 'Markersize', 4, 'MarkerFaceColor', 'black', 'Color', 'black');
    err.CapSize = 10;
    legends(4) = err;
    legend(legends, legendNames, 'Interpreter', 'latex', 'FontSize', 15, 'location', 'Southwest');
    title(['$\mathrm{Preamplifier \ Gain} \left( \mathrm{Gauge:} \ ',sprintf('%.1f', Gain),'\right)$'], 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('$\mathrm{Frequency[Hz]}$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$\mathrm{Gain[dB]}$', 'Interpreter', 'latex', 'fontsize', 20);
    hold off;
end