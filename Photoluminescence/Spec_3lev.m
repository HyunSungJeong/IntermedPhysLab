function I = Spec_3lev(I0, Linewidth, Ecen, Energy)
    % <Description>
    % Calculates normalized photoluminescence intensity at energies specified by 'Energy', 
    % according to the 3-level effective model
    %
    % <Input>
    % I0 : [numeric] The overall factor in the photoluminescence amplitude
    %
    % Linewidth : [numeric] Linewidth of the excited state
    % 
    % Ecen : [numeric] The energy difference between the ground state and the
    %                   first excited state |E1> in meV
    %
    % Energy : [numeric vector] The single-photon energies at which the
    %                           photoluminescence spectrum will be calculated
    %
    % <Output>
    % I : [numeric] Photoluminescence spectrum at energies specified by 'Energy'

    %% Parse input

    if ~isnumeric(I0) || ~isscalar(I0)
        error('ERR: ''I0'' must be a nonegative scalar');
    elseif I0 < 0
        error('ERR: ''I0'' must be nonnegative');
    end

    if ~isnumeric(Linewidth) || ~isscalar(Linewidth)
        error('ERR: ''I0'' must be a nonegative scalar');
    elseif Linewidth < 0
        error('ERR: ''I0'' must be nonnegative');
    end

    if ~isnumeric(Ecen) || ~isscalar(Ecen)
        error('ERR: ''E1'' must be a nonegative scalar');
    elseif Ecen < 0
        error('ERR: ''E1'' must be nonnegative');
    end

    if ~isnumeric(Energy)
        error('ERR: ''Energy'' must be a vector containing positive numbers');
    end

    %% Compute photoluminescence spectrum

    % physical constants
    el = 1.602 * 1e-19;     % elementary charge
    kb = 1.38 * 1e-23;      % Boltzmann constant
    h = 6.626 * 1e-34;      % Planck constant
    e_kb = 1e-3 * el / kb;

    % Intensity
    I = zeros(numel(Energy), 1);

    for it = 1:numel(Energy)
        I(it) = I0 * Lorenztian((Energy(it)-Ecen), Linewidth);
    end

    %% Define Lorentzain function
    function Lor = Lorenztian(X, Gamma)
        Lor = Gamma / (X^2 + (Gamma/2)^2) / 2;
    end

    function G = Gauss(X, sigma)
        G = exp(-X^2 / 2 / sigma^2) / sqrt(2*pi*sigma^2);
    end
end