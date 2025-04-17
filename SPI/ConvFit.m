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
    %            (I_0, lambda, thetaShift, Ishift) are fitting parameters, and there is no convolution
    %
    % 'FitLambda', ... : [numeric] Similar to 'NoConv1' option, but 'thetaShift' and 'Ishift' are fixed.
    %                            The inputs must be 'thetaShift' and 'Ishift', in this given order.
    %                            'thetaShift' must be in units of miliradians, and
    %                            'Ishift' must be either voltage(laser experiments) in mV units or # of photons(single-photon experiments).
    %                            (I_0, lambda) are fitting parameters.
    %
    % 'FitIntensity', ... : [numeric] Similar to 'NoConv2' option, but 'lambda' is also fixed.
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
    % 'Gaussian' : When used, the convolution of the chosen fitting model and the Gaussian function is used as the fitting model.
    %              The Gaussian function is defined as
    %                   Gauss(x) = (1/sqrt(2*pi*sigma^2))*exp(-x^2/(2*sigma^2))
    %              and sigma is a fitting parameter
    %              (Default: not used)
    %
    % <Output>
    % The fitting parameters are given as output.
    % The output parameters are naturally dependent on the fitting option used.
    % 1. 'NoConv' : [gamma, I_0, a, b, c, d]
    % 2. 'FitLorentz' : [gamma, I_0]
    % 3. 'FitSingle' : [gamma, I_0, a, b]
    % 4. 'FitDouble' : [gamma, I_0, b, c]
    % 5. 'FitAll' : [gamma, a, b, c]
    
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
    convFunc = NaN;

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
                ConvFunc = 'Lorentz';
                varargin(1) = [];

            case 'Gaussian'
                UseConv = true;
                ConvFunc = 'Gaussian';
                varargin(1) = [];

            otherwise
                if ~ischar(varargin{1})
                    error('ERR: Unknown input');
                else
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                end
        end
    end

    if isnan(fitModel)
        error('ERR: You must choose a fitting model');
    end
    
    %% Fit Data
    
    switch fitModel
        case 'FitAll'
            Params0 = [10000, 600, 0, 0];
            FitParams = lsqcurvefit(@NoConvModel, Params0, theta, I);
    
            if nargout ~= 4
                error('ERR: When ''FitAll'' fitting model is used without convolution, the number of output data must be 4');
            else
                varargout = num2cell(FitParams,1);
            end
    
        case 'FitLambda'

    
        case 'FitIntensity'
    end
    
    
    %% Define model functions

    % Lorentzian function
    function L = Lor(gamma, X)
        L = 1/(pi*gamma*(1 + X^2/gamma^2));
    end

    % \beta(theta) = 0.5*(2*pi/lambda)*slitWidth*sin(theta)
    function beta = Beta(lambda, slitWidth, theta)
        beta = 1e3*(pi/lambda)*slitWidth*sin(1e-3*theta);
    end

    % \gamma(theta) = 0.5*(2*pi/lambda)*slitDist*sin(theat)
    function gamma = Gamma(lambda, slitDist, theta)
        gamma = 1e3*(pi/lambda)*slitDist*sin(1e-3*theta);
    end

    % 'NoConv' model
    function I = NoConvModel(Params, theta)
        % Params = [I_0, lambda, thetaShift, Ishift]

        I = zeros(numel(theta),1);
        for it = 1:numel(I)
            if theta(it) ~= Params(3)
                I(it) = Params(1) * cos( Gamma(Params(2), slitDist, theta(it)-Params(3)) )^2;
                I(it) = I(it) * sin( Beta(Params(2), slitWidth, theta(it)-Params(3)) )^2 / Beta(Params(2), slitWidth, theta(it)-Params(3))^2;
                I(it) = I(it) + Params(4);
            else
                I(it) = Params(1) * cos( Gamma(Params(2), Params(3), theta(it)) )^2 + Params(4);
            end
        end
    end

    % 'FitLorentz' model
    function I = FitLorentzModel(Params, theta)

        I = zeros(numel(theta),1);
        Mono = @(x) NoConvModel([Params(1), a, b, c, 0], x); 
        for it = 1:numel(I)
            I(it) = Params(1)*integral();
        end
    end

end