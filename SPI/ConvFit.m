function varargout = ConvFit(theta, I, slitDist, slitWidth, varargin)
    % <Description>
    % Performs convolution fit of the input double-slit interference data
    %
    % <Input>
    % theta : [numeric vector] A vector containing angles from the screen center at which intensity(or number of photons) are measured.
    %                          'theta' must be in units of miliradians
    %
    % I : [numeric vector] A vector containing intensity or number of photons at angles given by the input variable 'theta'.
    %                      'I' must be either output voltage(laser experiments) in mV units or # of photons(single-photon experiments).
    %                      'theta' and 'I' must have same lengths.
    % slitDist : [numeric] Distance between slits in micrometets
    % slitWidth : [numeric] Slit width in micrometers
    %
    % <Options>
    %
    % 1. Fitting models: one of the three fitting models must be chosen
    %
    % 'FitAll' : If used, the data is fitted to the model
    %               I(theta) = I_0 * [sin( \beta(theta - thetaShift) ) / \beta(theta - thetaShift) ]^2 * cos^2(\gamma(theta - thetaShift)) + Ishift
    %            Where \beta(theta) = 0.5*(2*pi/lambda)*slitWidth*sin(theta)
    %                  \gamma(theta) = 0.5*(2*pi/lambda)*slitDist*sin(theat)
    %            (lambda, I_0, thetaShift, Ishift) are fitting parameters, and there is no convolution
    %
    % 'FitLambda', ... : [numeric] Similar to 'FitAll' option, but 'thetaShift' and 'Ishift' are fixed.
    %                            The inputs must be 'thetaShift' and 'Ishift', in this given order.
    %                            'thetaShift' must be in units of miliradians, and
    %                            'Ishift' must be either voltage(laser experiments) in mV units or # of photons(single-photon experiments).
    %                            (lambda, I_0) are fitting parameters.
    %
    % 'FitIntensity', ... : [numeric] Similar to 'FitLambda' option, but 'lambda' is also fixed.
    %                            The inputs must be 'lambda', 'thetaShift', and 'Ishift', in this given order.
    %                            'lambda' must be in units of nanometers and 'thetaShift' must be in units of miliradians.
    %                            'Ishift' must be # of photons(single-photon experiment) or voltage in mV(laser experment) units.
    %                            I_0 is the only fitting parameter.
    %
    % 2. Convolution fitting options
    %  Default : convolution fitting is not used
    %
    % 'Lorentz' : When used, the convolution of the chosen fitting model and the Lorentzian function is used as the fitting model.
    %             The Lorentzian function is defined as
    %                   Lor(x) = 1/(pi*gamma*(1 + x^2/gamma^2))
    %             and gamma is a fitting parameter
    %             (Default: not used)
    %
    % 'Gauss' : When used, the convolution of the chosen fitting model and the Gaussian function is used as the fitting model.
    %              The Gaussian function is defined as
    %                   Gauss(x) = (1/sqrt(2*pi*sigma^2))*exp(-x^2/(2*sigma^2))
    %              and sigma is a fitting parameter
    %              (Default: not used)
    %
    % 3. Other options
    %
    % 'InitParams', .. : [numeric vector] The initial value of the fitting parameters to be fed into 'lsqcurvefit'
    %
    % '-v' : When used, the data is plotted together with the fitted curve
    %           (Default: not used)
    %
    % <Output>
    % The fitting parameters are given as output.
    % The output parameters are naturally dependent on the fitting model used.
    % 1. 'FitAll' : [lambda, I_0, thetaShift, Ishift]
    % 2. 'FitLambda' : [lambda, I_0]
    % 3. 'FitIntensity' : [I_0]
    % 
    % If one of the convolution fitting options is used,
    % gamma('Lorentz' convolution) or sigma('Gauss' convolution) is appended to the first element of the output
    % e.g)
    % 'FitAll' and 'Lorentz' --> output is [gamma, lambda, I_0, thetaShift, Ishift]
    
    %% Parse inputs
    if ~isnumeric(theta)
        error('ERR: ''theta'' must be a vector of doubles');
    end
    
    if ~isnumeric(I)
        error('ERR: ''I'' must be a vector of doubles');
    end
    
    if ~isequal(numel(theta), numel(I))
        error('ERR: ''theta'' and ''I'' must be vectors of the same length');
    end

    if ~isnumeric(slitDist)
        error('ERR: ''slitDist'' must be a positive real number');
    elseif slitDist < 0
        error('ERR: ''slitDist'' must be positive');
    end

    if ~isnumeric(slitWidth)
        error('ERR: ''slitWidth'' must be a positive real number');
    elseif slitWidth < 0
        error('ERR: ''slitWidth'' must be positive');
    end
    
    % make sure fitting data are column vectors
    theta = theta(:);
    I = I(:);
    
    %% Parse options

    % Default options
    fitModel = NaN;
    ConvOp = NaN;
    UseConv = false;
    PlotFit = false;
    InitParamsGiven = false;

    % parse input options
    while ~isempty(varargin)
        switch varargin{1}
            case 'FitAll'
                fitModel = 'FitAll';
                varargin(1) = [];
    
            case 'FitLambda'
                fitModel = 'FitLambda';
    
                if ~isnumeric(varargin{2})
                    error('ERR: ''thetaShift'' must be a real number');
                elseif ~isnumeric(varargin{3})
                    error('ERR: ''Ishift'' must be a real number');
                else
                    thetaShift = varargin{2};
                    Ishift = varargin{3};
                end
                varargin(1:3) = [];
    
            case 'FitIntensity'
                fitModel = 'FitIntensity';
                
                if ~isnumeric(varargin{2})
                    error('ERR: ''lambda'' must be a positive real number');
                elseif varargin{2} < 0
                    error('ERR: ''lambda'' must be positive');
                elseif ~isnumeric(varargin{3})
                    error('ERR: ''thetaShift'' must be a real number');
                elseif ~isnumeric(varargin{4})
                    error('ERR: ''Ishift'' must be a real number');
                else
                    lambda = varargin{2};
                    thetaShift = varargin{3};
                    Ishift = varargin{4};
                end
                varargin(1:4) = [];
    
            case 'Lorentz'
                UseConv = true;
                ConvOp = 'Lorentz';
                varargin(1) = [];

            case 'Gauss'
                UseConv = true;
                ConvOp = 'Gauss';
                varargin(1) = [];

            case 'InitParams'

                if ~isnumeric(varargin{2})
                    error('ERR: Initial value of fitting parameters must be real numbers');
                else
                    InitParamsGiven = true;
                    Params0 = varargin{2};
                    varargin(1:2) = [];
                end

            case '-v'
                PlotFit = true;
                varargin(1) = [];

            otherwise
                if ~ischar(varargin{1})
                    error('ERR: Unknown input');
                else
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                end
        end
    end

    if InitParamsGiven
        if UseConv
            switch fitModel
                case 'FitAll'
                    if numel(Params0) ~= 5
                        error('ERR: The initial values of fitting parameters must be 5 in this fitting option');
                    end
    
                case 'FitLambda'
                    if numel(Params0) ~= 3
                        error('ERR: The initial values of fitting parameters must be 3 in this fitting option');
                    end
    
                case 'FitIntensity'
                    if numel(Params0) ~= 2
                        error('ERR: The initial values of fitting parameters must be 2 in this fitting option');
                    end
            end % switch-case
    
        else
            switch fitModel
                case 'FitAll'
                    if numel(Params0) ~= 4
                        error('ERR: The initial values of fitting parameters must be 4 in this fitting option');
                    end
    
                case 'FitLambda'
                    if numel(Params0) ~= 2
                        error('ERR: The initial values of fitting parameters must be 2 in this fitting option');
                    end
    
                case 'FitIntensity'
                    if numel(Params0) ~= 1
                        error('ERR: The initial values of fitting parameters must be 1 in this fitting option');
                    end
            end % switch-case
    
        end
    end

    if isnan(fitModel)
        error('ERR: You must choose a fitting model');
    end
    
    %% Fit Data
    
    switch fitModel
        case 'FitAll'       % If 'FitAll' model is used

            if isequal(ConvOp, 'Lorentz')       % convolution with Lorentzian
                Model = getConvModel(@FitAllModel, @Lor);
            elseif isequal(ConvOp, 'Gauss')     % convolution with Gaussian function
                Model = getConvModel(@FitAllModel, @Gauss);
            else        % without convolution
                Model = @FitAllModel;
            end

            if UseConv      % with convolution

                LowerBound = [0, 300, 0, -Inf, -Inf];
                UpperBound = [1e-2, 800, Inf(1,3)];

                if InitParamsGiven
                    for it = 1:numel(Params0)
                        if Params0(it) < LowerBound(it) || Params0(it) > UpperBound(it)
                            error('ERR: Initial value of fitting parameters is unphysical');
                        end
                    end
                else
                    Params0 = [1, 700, 1e4, 0, 0];
                end

                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);

                if nargout ~= 5
                    error('ERR: When ''FitAll'' fitting model is used with convolution, the number of output data must be 5');
                else
                    varargout = num2cell(FitParams,1);
                end

            else        % without convolution

                LowerBound = [300, 0, -Inf, -Inf];
                UpperBound = [800, Inf(1,3)];

                if InitParamsGiven
                    for it = 1:numel(Params0)
                        if Params0(it) < LowerBound(it) || Params0(it) > UpperBound(it)
                            error('ERR: Initial value of fitting parameters is unphysical');
                        end
                    end
                else
                    Params0 = [700, 1e4, 0, 0];
                end

                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);
        
                if nargout ~= 4
                    error('ERR: When ''FitAll'' fitting model is used without convolution, the number of output data must be 4');
                else
                    varargout = num2cell(FitParams,1);
                end
            end

    
        case 'FitLambda'

            FitLambdaModel = @(Params,theta) FitAllModel([Params(1:2),thetaShift,Ishift], theta);

            if isequal(ConvOp, 'Lorentz')       % convolution with Lorentzian
                Model = getConvModel(FitLambdaModel, @Lor);
            elseif isequal(ConvOp, 'Gauss')     % convolution with Gaussian function
                Model = getConvModel(FitLambdaModel, @Gauss);
            else        % without convolution
                Model = FitLambdaModel;
            end

            if UseConv      % with convolution

                LowerBound = [1e-8, 300, 2000];
                UpperBound = [1e-2, 800, 1e4];

                if InitParamsGiven
                    for it = 1:numel(Params0)
                        if Params0(it) < LowerBound(it) || Params0(it) > UpperBound(it)
                            error('ERR: Initial value of fitting parameters is unphysical');
                        end
                    end
                else
                    Params0 = [1, 700, 1e4];
                end
                
                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);
    
                if nargout ~= 3
                    error('ERR: When ''FitLambda'' fitting model is used with convolution, the number of output data must be 3');
                else
                    varargout = num2cell(FitParams,1);
                end
                

            else        % without convolution

                LowerBound = [300, 0];
                UpperBound = [800, Inf];

                if InitParamsGiven
                    for it = 1:numel(Params0)
                        if Params0(it) < LowerBound(it) || Params0(it) > UpperBound(it)
                            error('ERR: Initial value of fitting parameters is unphysical');
                        end
                    end
                else
                    Params0 = [700, 1e4];
                end

                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);
        
                if nargout ~= 2
                    error('ERR: When ''FitLambda'' fitting model is used without convolution, the number of output data must be 2');
                else
                    varargout = num2cell(FitParams,1);
                end

            end

        case 'FitIntensity'

            FitLambdaModel = @(Params,theta) FitAllModel([Params(1:2),thetaShift,Ishift], theta);
            FitIntensityModel = @(Params,theta) FitAllModel([lambda,Params(1),thetaShift,Ishift], theta);

            % ***we must use FitLambdaModel when convultion is used, since we have to perform k-integration
            if isequal(ConvOp, 'Lorentz')       % convolution with Lorentzian
                Model = getConvModel(FitLambdaModel, @Lor, 'FixLambda', lambda);
            elseif isequal(ConvOp, 'Gauss')     % convolution with Gaussian function
                Model = getConvModel(FitLambdaModel, @Gauss, 'FixLambda', lambda);
            else        % without convolution
                Model = FitIntensityModel;
            end

            if UseConv      % with convolution

                LowerBound = [0, 0];
                UpperBound = [1e-2, Inf];

                if InitParamsGiven
                    for it = 1:numel(Params0)
                        if Params0(it) < LowerBound(it) || Params0(it) > UpperBound(it)
                            error('ERR: Initial value of fitting parameters is unphysical');
                        end
                    end
                else
                    Params0 = [1e-1, 1.9*1e4];
                end

                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);
    
                if nargout ~= 2
                    error('ERR: When ''FitIntensity'' fitting model is used with convolution, the number of output data must be 2');
                else
                    varargout = num2cell(FitParams,1);
                end
                

            else        % without convolution

                LowerBound = 0;
                UpperBound = Inf;

                if InitParamsGiven
                    if Params0 < LowerBound || Params0 > UpperBound
                        error('ERR: Initial value of fitting parameter is unphysical');
                    end
                else
                    Params0 = 700;
                end

                FitParams = lsqcurvefit(Model, Params0, theta, I, LowerBound, UpperBound);
        
                if nargout ~= 1
                    error('ERR: When ''FitIntensity'' fitting model is used without convolution, the number of output data must be 1');
                else
                    varargout = num2cell(FitParams,1);
                end

            end

    end % switch-case

    if PlotFit      % if '-v' option is used
        %figure;
        %hold on;
        %plot(theta, I,'.','Color','black');
        thetaFit = linspace(theta(1), theta(end), 1e3);
        plot(thetaFit, Model(FitParams,thetaFit), 'Color', [.85 .325 .098]);
        %hold off;
    end
    
    
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