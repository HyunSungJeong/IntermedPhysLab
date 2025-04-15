function varargout = ConvFit(theta, I, varargin)
    % <Description>
    % Performs convolution fit of the input double-slit interference data
    %
    % <Input>
    % theta : [numeric vector] A vector containing angles from the screen center
    %                       at which intensity(or number of photons) are measured.
    % I : [numeric vector] A vector containing intensity or number of photons
    %                   at angles given by the input variable 'theta'.
    %                   'theta' and 'I' must have same lengths.
    %
    % <Options>
    % ***One option must be used, and the options are (obviously) not compatible
    % 'NoConv' : If used, the data is fitted to the model
    %               I = I_0 * [sin(a(x-b)) / (a(x-b))]^2 * cos^2[c(x-b)] + d
    %            Here, (I_0, a, b, c, d) are fitting parameters, and there is no convolution
    %
    % **The below fitting options are convolution fits, with different fitting parameters.
    % **The monocromatic fitting parameters (I_0, a, b, c, d) are defined as:
    %           I = I_0 * [sin(a(x-b)) / (a(x-b))]^2 * cos^2[c(x-b)] + d
    % **The width of the lorentzian, gamma, is defined as:
    %           Lor(lambda) = 1/(pi*gamma*(1 + (lambda - lambda_0)^2/gamma^2 ))
    %
    % 'FitLorentz', a, b, c : If used, gamma and I_0 are the fitting parameter.
    %                       a, b, c must be given as input after 'FitLorentz', and d is assumed to be zero.
    % 'FitSingle', c : If used, gamma and the single slit parameters (I__0, a, b) are the fitting parameters
    %               c must be given as input after 'FitSingle', and d is assumed to be zero.
    % 'FitDouble', a : If used, gamma and (I_0, b, c) are the fitting parameters.
    %               a must be given as input after 'FitDouble', and d is assumed to be zero.
    % 'FitAll' : If used, (gamma, a, b, c) are fitting parameters.
    %           d is assumed to be zero.
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
    
    theta = theta(:);
    I = I(:);
    
    %% Parse options
    while ~isempty(varargin)
        switch varargin{1}
            case 'NoConv'
                fitOption = 'NoConv';
                varargin(1) = [];
    
            case 'FitLorentz'
                fitOption = 'FitLorentz';
    
                if ~isnumeric(varargin{2})
                    error('ERR: ''a'' must be a real number');
                elseif ~isnumeric(varargin{3})
                    error('ERR: ''b'' must be a real number');
                elseif ~isnumeric(varargin{4})
                    error('ERR: ''c'' must be a real number');
                else
                    a = varargin{2};
                    b = varargin{3};
                    c = varargin{4};
                end
                varargin(1:4) = [];
    
            case 'FitSingle'
                fitOption = 'FitSingle';
                
                if ~isnumeric(varargin{2})
                    error('ERR: ''c'' must be a real number');
                else
                    c = varargin{2};
                end
                varargin(1:2) = [];
    
            case 'FitDouble'
                fitOption = 'FitDouble';
                
                if ~isnumeric(varargin{2})
                    error('ERR: ''a'' must be a real number');
                else
                    a = varargin{2};
                end
                varargin(1:2) = [];
    
            case 'FitAll'
                fitOption = 'FitAll';
                varargin(1) = [];
    
            otherwise
                if ~ischar(varargin{1})
                    error('ERR: Unknown input');
                else
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                end
        end
    end
    
    %% Fit Data
    
    switch fitOption
        case 'NoConv'
            Params0 = [1, 1, 0, 3, 0];
            FitParams = lsqcurvefit(@NoConvModel, Params0, theta, I);
    
            if nargout ~= 5
                error('ERR: When ''NoConv'' fitting option is used, the number of output data must be 5');
            else
                varargout = num2cell(FitParams,1);
            end
    
        case 'FitLorentz'
    
        case 'FitSingle'
    
        case 'FitDouble'
    
        case 'FitAll'
    end
    
    
    %% Model functions

    % 'NoConv' model
    function I = NoConvModel(Params, theta)

        I = zeros(numel(theta),1);
        for it = 1:numel(I)
            if theta(it) ~= Params(3)
                I(it) = Params(1)*( cos(Params(4)*(theta(it) - Params(3)))^2 );
                I(it) = I(it)*( sin(Params(2)*(theta(it) - Params(3))) / (Params(2)*(theta(it) - Params(3))) )^2;
                I(it) = I(it) + Params(5);
            else
                I(it) = Params(1)*cos(Params(4)*(theta(it) - Params(3)))^2 + Params(5);
            end
        end
    end

    % 'FitLorentz' model
    function I = FitLorentzModel(Params, theta)

        I = zeros(numel(theta),1);
        for it = 1:numel(I)
            
        end

    end

end