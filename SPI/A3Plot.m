% === File paths ===
path{1} = [fileparts(mfilename('fullpath')), filesep, 'spidata', filesep, 'D-Dark Data.txt'];
path{2} = [fileparts(mfilename('fullpath')), filesep, 'spidata', filesep, 'D-Small Light.txt'];
path{3} = [fileparts(mfilename('fullpath')), filesep, 'spidata', filesep, 'E1.txt'];
path{4} = [fileparts(mfilename('fullpath')), filesep, 'spidata', filesep, 'E2.txt'];

% === Initialize ===
position = cell(1,4);
I = cell(1,4);
STD = cell(1,4);

% === Load data ===
for i = 1:4
    data = importdata(path{i});
    if isstruct(data)
        data = data.data;
    end
    position{i} = data(:,1);
    I{i} = data(:,2);
    STD{i} = data(:,3);
end

position{3} = position{3} * 1000;
% === Define Seoul Metro Line 1–4 Colors ===
color{1} = [0.000, 0.322, 0.643]; % Line 1 (Dark Blue)
color{2} = [0.000, 0.659, 0.302]; % Line 2 (Green)
color{3} = [0.937, 0.486, 0.110]; % Line 3 (Orange)
color{4} = [0.000, 0.627, 0.914]; % Line 4 (Sky Blue)

% === Linear fit for Dark ===
p_dark = polyfit(position{1}, I{1}, 1);
yfit_dark = polyval(p_dark, position{1});
Rsq_dark = 1 - sum((I{1} - yfit_dark).^2) / sum((I{1} - mean(I{1})).^2);

% === Linear fit for Small Light ===
p_light = polyfit(position{2}, I{2}, 1);
yfit_light = polyval(p_light, position{2});
Rsq_light = 1 - sum((I{2} - yfit_light).^2) / sum((I{2} - mean(I{2})).^2);

% === Quadratic fit around peak of E1 ===
[~, peakIdx] = max(I{3});
range = max(1, peakIdx-5) : min(length(position{3}), peakIdx+5);
p_E1 = polyfit(position{3}(range), I{3}(range), 2);
xfit_E1 = position{3}(range);
yfit_E1 = polyval(p_E1, xfit_E1);
Rsq_E1 = 1 - sum((I{3}(range) - yfit_E1).^2) / sum((I{3}(range) - mean(I{3}(range))).^2);

% === Exponential fit for E2 ===
valid = I{4} > 0;
x_E2 = position{4}(valid);
logI = log(I{4}(valid));
p_E2 = polyfit(x_E2, logI, 1);
A = exp(p_E2(2));
B = p_E2(1);
xfit_E2 = linspace(min(position{4}), max(position{4}), 300); % Smooth x-grid
yfit_E2 = A * exp(B * xfit_E2);
Rsq_E2 = 1 - sum((logI - polyval(p_E2, x_E2)).^2) / sum((logI - mean(logI)).^2);

% === Labels ===
label = {'Operating Voltage - Dark Noise', 'Operating Voltage - Bulb Signal (Light Level 2)', 'Detector Position - Bulb Signal (Light Level 3)', 'Light Level - PCIT Gain'};

% === Plot each dataset separately ===
for i = 1:4
    figure;
    hold on; box on; % <-- No grid!

    % --- Plot error bars with smaller markers ---
    errorbar(position{i}, I{i}, STD{i}, 'o', ...
        'Color', 'k', 'MarkerFaceColor', 'k', ...
        'MarkerSize', 4, 'LineWidth', 1.2, ...
        'DisplayName', 'PCIT Gain');

    % --- Plot fitted curve ---
    switch i
        case 1 % Dark
            plot(position{1}, yfit_dark, '-', 'Color', color{i}, 'LineWidth', 2, 'DisplayName', 'Linear Fit');
            eqn = {sprintf('y = %.3fx %.3f', p_dark(1), p_dark(2)), ...
                   sprintf('R^2 = %.4f', Rsq_dark),sprintf('Threshold V = %.1f',(-p_dark(2)/p_dark(1)))};
        case 2 % Light
            plot(position{2}, yfit_light, '-', 'Color', color{i}, 'LineWidth', 2, 'DisplayName', 'Linear Fit');
            eqn = {sprintf('y = %.3fx %.3f', p_light(1), p_light(2)), ...
                   sprintf('R^2 = %.4f', Rsq_light),sprintf('Threshold V = %.1f',(-p_light(2)/p_light(1)))};
        case 3 % E1
            plot(xfit_E1, yfit_E1, '-', 'Color', color{i}, 'LineWidth', 2, 'DisplayName', 'Quadratic Fit (Peak)');
            eqn = {sprintf('peak at x = %.3f', - p_E1(2) / (2*p_E1(1))), ...
                   sprintf('R^2 = %.4f', Rsq_E1)};
        case 4 % E2
            plot(xfit_E2, yfit_E2, '-', 'Color', color{i}, 'LineWidth', 2, 'DisplayName', 'Exponential Fit');
            eqn = {sprintf('y = %.3fe^{%.3fx}', A, B), ...
                   sprintf('R^2 = %.4f', Rsq_E2)};
    end

    % --- Adjust axes ---
    dx = 0.05 * (max(position{i}) - min(position{i}));
    xlim([min(position{i}) - dx, max(position{i}) + dx]);

    if i <= 2
        % For Dark and Light (former 2 datasets): Fixed y-limits
        ylim([-2, 2]);
        xlabel('PMT Operating Voltage (V)','FontSize',25);
        ylabel('log(PCIT Gain/s)','FontSize',25);
        if i == 1
            xline(-p_dark(2)/p_dark(1),'-r','LineWidth',1.0,'DisplayName','Threshold V');
        end
        if i==2
            xline(-p_light(2)/p_light(1),'-r','LineWidth',1.0,'DisplayName','Threshold V');
        end
    else
        % For E1 and E2: Dynamic y-limits with margin
        dy = 0.1 * (max(I{i}) - min(I{i}));
        ylim([min(I{i}) - dy, max(I{i}) + dy]);
        ylabel('PCIT Gain/s','FontSize',25);
        if i==3
            xlabel('Position of the Detector Slit (mm)','FontSize',25);
        end
        if i==4
            xlabel('Light Level of the Bulb','FontSize',25);
        end
    end

    % --- Draw y=0 guideline ---
    yline(0, '--k', 'LineWidth', 1.0, 'DisplayName','Gain per Second = 1'); % Dashed black line at y=0

    % --- Annotate fitting equation ---
    text(mean(xlim), 0.8*max(ylim), eqn, ...
         'Color', color{i}, 'FontSize', 20, 'HorizontalAlignment', 'center');
    title(sprintf('%s', label{i}),"FontSize",30);
    legend('Location', 'best');
    set(gca, 'FontSize', 25);
    hold off;
end

