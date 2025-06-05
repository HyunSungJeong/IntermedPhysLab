function getPhaseShift_fixedPhase()
    Delta = zeros(1, 16);
    PhaseStd = zeros(1, 16);
    Freq = [10, 75, 100*2.^(0:13)];

    for it = 3:16
        phaseShift = 30;  % Fixed gauge shift
        dominant_freq = Freq(it);

        % --------- Step 1: Read the CSV and find the header row ---------
        if it < 3
            path = [fileparts(mfilename('fullpath')), filesep, 'data', filesep, 'A2_fixedPhase', filesep, sprintf('%d', Freq(it)), 'Hz_PhaseShift30.xlsx'];
        else
            path = [fileparts(mfilename('fullpath')), filesep, 'data', filesep, 'A2_fixedPhase', filesep, sprintf('%d', Freq(it)), 'Hz_PhaseShift30.csv'];
        end

        raw = readcell(path);
        startRow = find(strcmpi(raw(:,1), 'TIME'), 1, 'first') + 1;
        data = readmatrix(path, 'NumHeaderLines', startRow - 1);

        if it < 3
            time = data(:,1);
            sig1 = data(:,2);
            sig2 = data(:,4);
        else
            time = data(:,1);
            sig1 = data(:,2);
            sig2 = data(:,3);
        end

        fs = 1 / mean(diff(time));  % Sampling rate

        % Calculate samples per cycle and number of full cycles
        samples_per_cycle = round(fs / dominant_freq);
        num_cycles = floor(length(sig1) / samples_per_cycle);

        % Trim to full cycles
        sig1 = sig1(1:num_cycles*samples_per_cycle);
        sig2 = sig2(1:num_cycles*samples_per_cycle);

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

            % Cross-correlation
            [c_val, lags] = xcorr(seg2, seg1, 'coeff');
            [~, max_idx] = max(c_val);
            lag = lags(max_idx);

            delta_t = lag / fs;
            %phase_shifts(c) = mod(delta_t * dominant_freq * 360, 360);  % In degrees
            phase_shifts(c) = delta_t * dominant_freq * 360;
        end

        % Compute mean and std deviation
        mean_phase = mean(phase_shifts);
        std_phase = std(phase_shifts);

        % Store anomalous phase shift (in degrees)
        Delta(it) = mean_phase - phaseShift;
        PhaseStd(it) = std_phase;

        fprintf('File %02d: Phase = %.2f ± %.2f deg (Expected: %.1f, Frequency: %d Hz)\n', it, mean_phase, std_phase, phaseShift, Freq(it));
    end

    % Plot results
    figure;
    hold on;
    ax = gca;
    ax.XLabel.FontSize = 15;
    ax.YLabel.FontSize = 15;
    set(gca, 'XScale', 'log');
    err = errorbar(Freq(3:16), deg2rad(Delta(3:16)), deg2rad(PhaseStd(3:16)), 'o', 'Markersize', 4, 'MarkerFaceColor', 'black', 'Color', 'black');
    err.CapSize = 10;
    title('$\mathrm{Anomalous \ Phase \ Shift \ at \ Gauge \ 30}$', 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('$\mathrm{Frequency [Hz]}$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$\mathrm{Anomalous \ Phase \ Shift [rad]}$', 'Interpreter', 'latex', 'fontsize', 20);
    grid on;
    hold off;
end