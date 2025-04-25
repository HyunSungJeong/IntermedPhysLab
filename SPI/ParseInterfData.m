function ParseInterfData()

    path = [fileparts(mfilename('fullpath')),filesep,'data',filesep,'LASER_Double_14.txt'];
    slitDist = 356;
    slitWidth = 85;


    data = importdata(path);
    if isstruct(data)
        data = data.data;
    end

    theta = 1e3*data(:,1);
    I = data(:,2);

    %{
    Params0 = [600, 1500, 4, 0];
    [lambda, I_0, thetaShift, Ishift] = ConvFit(theta,I,slitDist,slitWidth,'FitAll','InitParams',Params0,'-v');
    disp([lambda, I_0, thetaShift, Ishift]);
    %}

    %{
    Params0 = [1e-5, 600, 1500, 4, 0];
    [gamma, lambda, I_0, thetaShift, Ishift] = ConvFit(theta,I,slitDist,slitWidth,'FitAll','Lorentz','InitParams',Params0,'-v');
    disp(['gamma = ',sprintf('%.15g',gamma)]);
    disp([lambda, I_0, thetaShift, Ishift]);
    %}

    %{}
    Params0 = [1e-5, 600, 1500];
    [gamma, lambda, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitLambda',4.4,0,'Lorentz','InitParams',Params0,'-v');
    disp(['gamma = ',sprintf('%.15g',gamma)]);
    disp([lambda, I_0]);
    %}

    %{
    Params0 = [600, 1500];
    [lambda, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitLambda',4.4,0,'InitParams',Params0,'-v');
    disp([lambda, I_0]);
    %}

    %{
    Params0 = [1e-5, 1500];
    [gamma, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitIntensity',670,4.4,0,'Lorentz','InitParams',Params0,'-v');
    disp(['gamma = ',sprintf('%.15g',gamma)]);
    disp(I_0);
    %}




    %% Define frequently used functions

    % \beta(theta) = 0.5*(2*pi/lambda)*slitWidth*sin(theta)
    function beta = Beta(lambda, slitWidth, theta)
        beta = 1e3*(pi/lambda)*slitWidth*sin(1e-3*theta);
    end

    % \gamma(theta) = 0.5*(2*pi/lambda)*slitDist*sin(theta)
    function gamma = Gamma(lambda, slitDist, theta)
        gamma = 1e3*(pi/lambda)*slitDist*sin(1e-3*theta);
    end

    % Lorentzian function
    function L = Lor(gamma, X)
        L = 1./(pi*gamma*(1 + X.^2/gamma^2));
    end

    % Gaussian function
    function G = Gauss(sigma, X)
        G = (1/sqrt(2*pi*sigma^2))*exp(-X.^2/(2*sigma^2));
    end

    %% Define monocromatic interference function

    % 'FitAll' model, without convolution
    function I = FitAllModel(Params, theta)
        % Params = [lambda, I_0, thetaShift, Ishift]

        I = zeros(numel(theta),1);
        for itx = 1:numel(I)
            if theta(itx) ~= Params(3)
                I(itx) = Params(2) * cos( Gamma(Params(1), slitDist, theta(itx)-Params(3)) )^2;
                I(itx) = I(itx) * sin( Beta(Params(1), slitWidth, theta(itx)-Params(3)) )^2 / Beta(Params(1), slitWidth, theta(itx)-Params(3))^2;
                I(itx) = I(itx) + Params(4);
            else
                I(itx) = Params(2) * cos( Gamma(Params(1), slitDist, theta(itx)-Params(3)) )^2 + Params(4);
            end
        end
    end


    %% Define convolution function

    function ConvModelHandle = getConvModel(ModelFunc, ConvKer, varargin)
        % <Description>
        % Generates function handle for the convolution of given monocromatic function and the convolution kernel
        %
        % <Input>
        % ModelFunc : function hande for the model before convolution
        %             The first element of the fitting parameter of ModelFunc, Params(1), must be lambda(wavelength of light in nanometers)
        % ConvKer : function handle for the convolution kernel
        %
        % <Option>
        % 'FixLambda', ... : [numeric] If the central value of lambda is specified by this options,
        %                              lambda is treated as a fixed value instead of a fitting parameter
        %                               (Default: not used. Central value of lambda is fitted from data)
        %
        % <Output>
        % Function handle of the model function obtained by convolution of ModelFunc and ConvFunc
        % Let Params_orig be the fitting parameters of ModelFunc, and let Params_conv be the fitting parameter of the convolution kernel.
        % Then, Params = [Params_orig, Params_conv] is the fitting parameter of the output function, ConvFunc.
        % Therefore, ConvFunc is a function of Params = [Params_orig, Params_conv] and theta.

        lambdaFixed = false;

        while ~isempty(varargin)
            switch varargin{1}
                case 'FixLambda'
                    if ~isnumeric(varargin{2})
                        error('ERR: lambda must be a positive real number');
                    elseif varargin{2} < 0
                        error('ERR: lambda must be positive');
                    else
                        lambda_C = varargin{2};
                        lambdaFixed = true;
                        varargin(1:2) = [];
                    end
                otherwise
                    if ~ischar(varargin{1})
                    error('ERR: Unknown input for getConvModel');
                    else
                        error(['ERR: Unknown option ''',varargin{1},''' for getConvModel']);
                    end
            end
        end

        function I = ConvModel(Params, theta)

            function Y = ModelFunc_k(k,theta)   % ModelFunc in k-space
                Y = zeros(size(k));
                for itx = 1:numel(Y)
                    if ~lambdaFixed
                        Y(itx) = ModelFunc([2*pi/k(itx), Params(3:end)], theta);
                    else
                        Y(itx) = ModelFunc([2*pi/k(itx), Params(2:end)], theta);
                    end
                end
            end

            I = zeros(numel(theta),1);
            for ity = 1:numel(I)

                % define convolution integrand
                if ~lambdaFixed     % If lambda is not fixed by input
                    ConvInt = @(k) ModelFunc_k(k+2*pi/Params(2), theta(ity)).*ConvKer(Params(1), k);        
                else                % If lambda is fixed by input
                    ConvInt = @(k) ModelFunc_k(k+2*pi/lambda_C, theta(ity)).*ConvKer(Params(1), k);
                end

                I(ity) = Integrate_symlog(ConvInt, -10*Params(1), 10*Params(1), 2e2*(16+log10(Params(1))) );     % convolution integration
            end
        end

        ConvModelHandle = @ConvModel;
    end

    %% Define symmetric-log integration function

    function Int = Integrate_symlog(f, min, max, N)
        % <Description>
        % numerically integrates f(x) from min to max(min < 0 < max) using a symmetric logarithmic grid around zero
        %
        % <Input>
        % f : [function handle] function handle of the integrand.
        %                       The integrand is assumed to be sharply peaked around zero.
        % min : [numeric] lower limit of the integral. Must be negative
        % max : [numeric] upper limit of the integral. Must be positive
        % N : [numeric] number of total logarithmic grid points.
        %
        % <Output>
        % Int : [numeric] approximate value of the integral \int_{min}^{max} f(x) dx
        
        if min >= 0 || max <= 0
            error('ERR: This function requires min < 0 and max > 0 to integrate around zero.');
        end
    
        if mod(N, 2) ~= 0
            N = N + 1; % Make N even
        end
        N_half = N / 2;
    
        xL = - logspace(log10(abs(min)), log10(1e-16), N_half); % create log grid for left(negative) side
        xR = logspace(log10(1e-16), log10(max), N_half);        % create log grid for right(positive) side
        x = [xL, xR];   % combine grid   
        y = f(x);       % evaluate function on grid
    
        Int = trapz(x, y);      % integrate using trapz
    end
end




function [k_dense, g_smooth] = smoothen_symlog_g(k, g_weights, N_dense)
    % SMOOTHEN_SYMLOG_G: Smoothens g(k) on a symmetric log scale
    %
    % Inputs:
    %   k          - original symmetric log grid (e.g., [-logspace(...), logspace(...)])
    %   g_weights  - recovered g(k) values on original k
    %   N_dense    - number of points per side (optional, default = 200)
    %
    % Outputs:
    %   k_dense    - finer symmetric log grid
    %   g_smooth   - smoothened g(k) values on finer grid

    if nargin < 3
        N_dense = 200;
    end

    % Build finer symmetric log grid
    k_min = min(abs(k(k ~= 0)));
    k_max = max(abs(k));
    k_dense_pos = logspace(log10(k_min), log10(k_max), N_dense);
    k_dense = [-fliplr(k_dense_pos), k_dense_pos];

    % Fit a smoothing spline (optional: control smoothness with p)
    % p = 0: interpolating spline, p = 1: least-squares linear fit
    p = 1e-4;  % You can adjust this between [0, 1]
    spline_fit = fit(k(:), g_weights(:), 'smoothingspline', 'SmoothingParam', p);

    % Evaluate on new grid
    g_smooth = feval(spline_fit, k_dense);
end
