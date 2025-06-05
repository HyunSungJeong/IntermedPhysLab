function getPhaseShift_fixedFreq(range)

    Delta = zeros(1, 24);
    PhaseStd = zeros(1, 24);

    % Parameters
    fs = 1e9;              % Sampling rate
    f = 1e5;               % Signal frequency (Hz)
    samples_per_cycle = round(fs / f);
    num_cycles = 10;       % Analyze 10 cycles per file

    for it = range
        phaseShift = 15*(it-1);

        % --------- Step 1: Read the CSV and find the header row ---------
        path = [fileparts(mfilename('fullpath')),filesep,'data',filesep,'A2_fixedFreq',filesep,'100kHz_',sprintf('%d',phaseShift),'Deg.csv'];

        raw = readcell(path);
        startRow = find(strcmpi(raw(:,1), 'TIME'), 1, 'first') + 1;
        data = readmatrix(path, 'NumHeaderLines', startRow - 1);

        % Columns: time | signal1 | signal2
        sig1 = data(:,2);
        sig2 = data(:,3);

        % Trim to exact number of samples
        N = samples_per_cycle * num_cycles;
        sig1 = sig1(1:N);
        sig2 = sig2(1:N);

        % Remove DC offset
        sig1 = sig1 - mean(sig1);
        sig2 = sig2 - mean(sig2);

        % Initialize per-cycle phase shift array
        phase_shifts = zeros(1, num_cycles);

        for c = 1:num_cycles
            idx_start = (c-1)*samples_per_cycle + 1;
            idx_end = c*samples_per_cycle;

            seg1 = sig1(idx_start:idx_end);
            seg2 = sig2(idx_start:idx_end);

            % Cross-correlation per cycle
            [c_val, lags] = xcorr(seg2, seg1, 'coeff');
            [~, max_idx] = max(c_val);
            lag = lags(max_idx);

            % Phase shift in degrees
            delta_t = lag / fs;
            phase_shifts(c) = mod(delta_t * f * 360, 360);
        end

        % Compute stats
        mean_phase = mean(phase_shifts);
        std_phase = std(phase_shifts);

        % Store anomalous phase shift relative to expected
        Delta(it) = deg2rad(mean_phase - phaseShift);
        PhaseStd(it) = deg2rad(std_phase);

        % Output
        fprintf('File %02d: Phase = %.2f ± %.2f deg (Expected: %d)\n', it, mean_phase, std_phase, phaseShift);
    end

    % Phase shift array for plotting
    phaseShift = 0:15:345;

    % Plot 1: Anomalous phase shift
    figure;
    hold on;
    ax = gca;
    ax.XLabel.FontSize = 15;
    ax.YLabel.FontSize = 15;
    err = errorbar(phaseShift(range), Delta(range), PhaseStd(range), 'o', 'Markersize', 4, 'MarkerFaceColor', 'black', 'Color', 'black');
    err.CapSize = 10;
    title('$\mathrm{Anomalous \ Phase \ Shift \ at \ 100 kHz}$', 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('$\mathrm{Gauge \ Phase \ Shift [rad]}$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$\mathrm{Anomalous \ Phase \ Shift [rad]}$', 'Interpreter', 'latex', 'fontsize', 20);
    hold off;

    % Plot 2: Measured phase shift
    figure;
    hold on;
    ax = gca;
    ax.XLabel.FontSize = 15;
    ax.YLabel.FontSize = 15;
    phaseShift_meas = deg2rad(phaseShift(range)) + Delta(range);
    err = errorbar(phaseShift(range), phaseShift_meas, PhaseStd(range), 'o', 'Markersize', 4, 'MarkerFaceColor', 'black', 'Color', 'black');
    err.CapSize = 10;
    title('$\mathrm{Phase \ Shift \ at \ 100 kHz}$', 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('$\mathrm{Gauge \ Phase \ Shift [rad]}$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$\mathrm{Measured \ Phase \ Shift [rad]}$', 'Interpreter', 'latex', 'fontsize', 20);
    hold off;
end
