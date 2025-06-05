clear;
freq_6 = {[0.1, 0.5, 1, 5, 10, 50, 100, 500];
            [0.1, 0.5, 1, 5, 10, 50, 100]};
Gain_6 = {[1.04717, 0.994382, 1.015694, 0.731771, 0.5, 0.157269, 0.103156, 0.086449];
        [1.037383, 0.990654, 0.878505, 0.364486, 0.226415, 0.103298, 0.093752]};
        
Gain_stdev_6 = {[0, 0, 0.005498, 0.004823, 0, 0.008747, 0.00988, 0.014019];
                [0, 0, 0, 0, 0, 0.001252, 0.00936]};

freq_12 = [0.5, 1, 5, 10, 50];
Gain_12 = {[0.994382, 1.004761, 0.542056, 0.252336, 0.093458];
            [1, 0.747664, 0.158879, 0.102804, 0.084685]};

%% plot Roll-off 6

figure;
hold on;
grid on;
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
ax.XAxis.MinorTick = 'on';
set(gca, 'XScale', 'log');

plot1 = errorbar(freq_6{1}, Gain_6{1}, Gain_stdev_6{1}, 'o', 'Markersize', 4, 'MarkerFaceColor', 'blue', 'Color', 'blue');
plot1.CapSize = 10;
plot2 = errorbar(freq_6{2}, Gain_6{2}, Gain_stdev_6{2}, 'o', 'Markersize', 4, 'MarkerFaceColor', 'red', 'Color', 'red');
plot2.CapSize = 10;
plot(freq_6{1}, Gain_6{1}, 'LineWidth', 0.7, 'Color', 'blue');
plot(freq_6{2}, Gain_6{2}, 'LineWidth', 0.7, 'Color', 'red');

legendNames = {'$\mathrm{time \ constant \ = \ 0.03 }$', '$\mathrm{time \ constant \ = \ 0.1 }$'};
legend([plot1, plot2], legendNames, 'Interpreter', 'latex', 'FontSize', 15, 'location', 'northeast');
title('$\mathrm{Roll-off: \ 6dB/oct}$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$\mathrm{Frequency[Hz]}$', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\mathrm{Gain[dB]}$', 'Interpreter', 'latex', 'fontsize', 20);
hold off;

%% plot Roll-off 12

figure;
hold on;
grid on;
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
ax.XAxis.MinorTick = 'on';
set(gca, 'XScale', 'log');

plot1 = plot(freq_12, Gain_12{1}, 'o', 'MarkerSize', 4, 'Color', 'blue');
plot2 = plot(freq_12, Gain_12{2}, 'o', 'MarkerSize', 4, 'Color', 'red');
plot(freq_12, Gain_12{1}, 'LineWidth', 0.7, 'Color', 'blue');
plot(freq_12, Gain_12{2}, 'LineWidth', 0.7, 'Color', 'red');

legendNames = {'$\mathrm{time \ constant \ = \ 0.03 }$', '$\mathrm{time \ constant \ = \ 0.1 }$'};
legend([plot1, plot2], legendNames, 'Interpreter', 'latex', 'FontSize', 15, 'location', 'northeast');
title('$\mathrm{Roll-off: \ 12dB/oct}$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$\mathrm{Frequency[Hz]}$', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\mathrm{Gain[dB]}$', 'Interpreter', 'latex', 'fontsize', 20);
hold off;