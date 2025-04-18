function varargout = ParseViscosData(N,D,T,varargin)
    % <Description>
    % Parses raw Tracker data from Brownian motion experiment and gives 
    % processed data such as x, y displacement, average of distance sqaured, etc.
    %
    % <Input>
    % N : [integer] Number of tracked particles data
    % D : [numeric] Diameter of microscopic particle used in experiment, in units of micrometer
    % T : [numeric] Temperature in Kelvins
    % 
    % <Options>
    % 'MaxNumData', ... : [integer] maximum number of sequential data for each particle
    %       (Default: 3000)
    % 'CorrectDrift' ... : If used, the average drift of the particles are subtracted from the original data
    % 
    %                      When used together with a input of 1xN cell array of two dimensional vectors,
    %                      These vectors are regarded drift velocities of each particle,
    %                      and are used to correct drift of individual particles
    %       (Default: not used)
    % 'PlotTimeRange' ... : [numeric] The time range to be used for fitting and plotting is given by [0, PlotTimeRange]
    %       (Default: 120)
    % '-v' : If used, the <r^2> vs time is plotted together wit a linear fit
    %       (Default: not used)
    % 'ShowText' : If used together with '-v' option,
    %           viscosity, R^2 value of linear fit, and the fitted equation is shown together with the plot
    %       (Default: not used)
    %
    % <Output>
    % There are three possible outputs:
    % 1. [Time,X,Y,Viscosity]
    % 2. [Time,X,Y,Viscosity,time,DistSqAvg,DistSqVar]
    % 3. [Time,X,Y,Viscosity,time,DistSqAvg,DistSqVar,DriftVel] : can only be used when 'CorrectDrift' option is on
    %
    % Time : [cell array of numeric vectors] Each cell element is a vector of frame times, for each particle. 
    % X : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the x-displacement of each particle 
    %           at times given by output variable 'time'
    % Y : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the y-displacement of each particle 
    %           at times given by output variable 'time'
    % Viscosity : [numeric] Calculated viscosity of the medium
    % time : [numeric vector] time for DistSqAvg and DistSqVar
    % DistSqAvg : [numeric vector] Numeric vector of the averages of displacemet squared at times give by output variable 'Time'
    % DistSqVar : [numeric vector] Numeric vector of the variance of displacement squared at times given by output variable 'Time'

    %% Parse inputs
    if ~isnumeric(N)
        error('ERR: ''N'' must be a positive integer');
    elseif N <= 0
        error('ERR: ''N'' must be a positive integer');
    end

    if ~isnumeric(D)
        error('ERR: ''D'' must be a positive real number');
    elseif D <= 0
        error('ERR: ''D'' must be a positive real number');
    end

    if ~isnumeric(T)
        error('ERR: ''T'' must be a positive real number');
    elseif T <= 0
        error('ERR: ''T'' must be a positive real number');
    end

    %% Parse options

    % Default values of options
    plotFig = false;
    CorrectDrift_default = false;
    CorrectDrift_custom = false;
    MaxNumData = 3000;
    PlotTimeRange = 120;
    ShowText = false;

    while ~isempty(varargin)
        switch varargin{1}
            case '-v'
                plotFig = true;
                varargin(1) = [];

            case 'CorrectDrift'
                if numel(varargin) > 1
                    if iscell(varargin{2})
                        CorrectDrift_custom = true;
                        DriftVel = varargin{2};
                        varargin(1:2) = [];
                    else
                        CorrectDrift_default = true;
                        varargin(1) = [];
                    end
                else
                    CorrectDrift_default = true;
                    varargin(1) = [];
                end

            case 'MaxNumData'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        MaxNumData = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''MaxNumData'' must be a positive integer');
                    end
                else
                    error('ERR: ''MaxNumData'' must be a positive integer');
                end

            case 'PlotTimeRange'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        PlotTimeRange = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''PlotTimeRange'' must be a positive number');
                    end
                else
                    error('ERR: ''PlotTimeRange'' must be a positive number');
                end

            case 'ShowText'
                ShowText = true;
                varargin(1) = [];

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end
    end

    if ~plotFig & ShowText
        disp('WRN: ''ShowText'' option is redundant without ''-v'' option');
    end


    %% Parse raw data

    time = linspace(0,(MaxNumData-1)/10,MaxNumData);
    kb = 1.38*1e-23;
    
    X = cell(1,N);
    Y = cell(1,N);
    DistSq = cell(1,N);
    DataIdx = cell(1,N);
    
    for it = 1:N
        % TODO
        % Change the path to the folder where data are stored
        data = importdata([fileparts(mfilename('fullpath')),filesep,'data',filesep,sprintf('%.15g',D),'um',filesep,'Viscosdata',num2str(it),'.txt']);   % load data
        
        if isstruct(data)
            data = data.data;
        end
        DataIdx{it} = nan(1,size(data,1));
        X{it} = nan(1,MaxNumData);
        Y{it} = nan(1,MaxNumData);
        DistSq{it} = nan(1,MaxNumData);
    
        skipIdx = [];
        for itD = 1:size(data,1)
            Time_now = data(itD,1);
            X_now = data(itD,2);
            Y_now = data(itD,3);
            
            if itD > 1
                X_prev = data(itD-1,2);
                Y_prev = data(itD-1,3);
                Time_prev = data(itD-1,1);
    
                if sqrt((X_now-X_prev)^2 + (Y_now-Y_prev)^2) > 1e-6     % if there is a sudden jump of coordinates
        
                    % calculate the amount of sudden jump
                    TimeJump = Time_now - Time_prev;
                    Xjump = X_now - X_prev;
                    Yjump = Y_now - Y_prev;
        
                    skipIdx = [skipIdx, itD];   % indices to be skipped when processing data into displacements

                    % subtract sudden jump from following data
                    data(itD:end,1) = data(itD:end,1)-  repmat(TimeJump,[size(data,1)-itD+1, 1]);
                    data(itD:end,2) = data(itD:end,2) - repmat(Xjump,[size(data,1)-itD+1, 1]);
                    data(itD:end,3) = data(itD:end,3) - repmat(Yjump,[size(data,1)-itD+1, 1]);
                end
            end
        end % itD
    
        for itD = 1:size(data,1)        % process data into displacements
            if ~ismember(itD,skipIdx)   % if there was no sudden jump in this frame
                DataIdx{it}(itD) = round(10*data(itD,1)) + 1;
                X_tmp = 1e6*(data(itD,2) - data(1,2));
                Y_tmp = 1e6*(data(itD,3) - data(1,3));
                X{it}(DataIdx{it}(itD)) = X_tmp;
                Y{it}(DataIdx{it}(itD)) = Y_tmp;
            end
        end
    
    end % it
    

    %% Correct overall drift of particles

    DriftAvg_X = 0;
    DriftAvg_Y = 0;
    Drift_num = 0;
    for it = 1:N
        if ~isnan(X{it}(round(10*PlotTimeRange)))
            DriftAvg_X = DriftAvg_X + X{it}(round(10*PlotTimeRange));
            DriftAvg_Y = DriftAvg_Y + Y{it}(round(10*PlotTimeRange));
            Drift_num = Drift_num + 1;
        else
            error(['ERR: some tracking data do not have data at time',sprintf('%.15g',PlotTimeRange),'. Change ''PlotTimeRange''']);
        end
    end
    DriftAvg_X = DriftAvg_X/Drift_num;
    DriftAvg_Y = DriftAvg_Y/Drift_num;
    if CorrectDrift_default
        DriftVel = [DriftAvg_X/PlotTimeRange, DriftAvg_Y/PlotTimeRange];
    elseif ~CorrectDrift_custom
        DriftVel = [0,0];
    end
    
    
    for it = 1:N
        if CorrectDrift_custom
            X{it} = X{it} - DriftVel{it}(1)*time;
            Y{it} = Y{it} - DriftVel{it}(2)*time;
            disp(DriftVel{it});
        else
            X{it} = X{it} - DriftVel(1)*time;
            Y{it} = Y{it} - DriftVel(2)*time;
        end
        DistSq{it} = X{it}.^2 + Y{it}.^2;
    end
    
    DistSqAvg = nan(1,MaxNumData);
    DistSqVar = nan(1,MaxNumData);
    for itD = 1:MaxNumData
        cnt = 0;
        for it = 1:N
            if ~isnan(DistSq{it}(itD))
                if isnan(DistSqAvg(itD)) && isnan(DistSqVar(itD))
                    DistSqAvg(itD) = DistSq{it}(itD);
                    DistSqVar(itD) = DistSq{it}(itD)^2;
                else
                    DistSqAvg(itD) = DistSqAvg(itD) + DistSq{it}(itD);
                    DistSqVar(itD) = DistSqVar(itD) + DistSq{it}(itD)^2;
                end
                cnt = cnt + 1;
            end
        end
        if cnt > 0
            DistSqAvg(itD) = DistSqAvg(itD)/cnt;
            DistSqVar(itD) = DistSqVar(itD)/cnt - DistSqAvg(itD)^2;
        end
    end
    
    %% Calculate viscosity from processed data

    % discard frames without data
    
    Time = cell(1,N);
    for it = 1:N
        Time{it} = time;
        EmptyIdx = isnan(X{it});
        Time{it}(EmptyIdx) = [];
        X{it}(EmptyIdx) = [];
        Y{it}(EmptyIdx) = [];
        DistSq{it}(EmptyIdx) = [];
    end

    PlotFrames = round(10*PlotTimeRange);   % Number of frames to be used for fitting and plotting
    if PlotFrames > numel(DistSqAvg)
        disp(numel(DistSqAvg))
        error('ERR: ''PlotTimeRange'' must be smaller than the time length of tracked data');
    end
    
    EmptyIdx = isnan(DistSqAvg);
    DistSqAvg(EmptyIdx) = [];
    DistSqVar(EmptyIdx) = [];
    time(EmptyIdx) = [];

    % linear fit <r^2> - time data and calculate viscosity
    mdl = fitlm(time(1:PlotFrames),DistSqAvg(1:PlotFrames),'Intercept',false);
    slope = mdl.Coefficients.Estimate;
    Rsq = mdl.Rsquared.Ordinary;
    coeff = [slope,0];

    Viscosity = (1e18)*4*kb*T/(3*pi*D*slope);
    Stdev = sqrt(DistSqVar);


    %% Plot figure if '-v' option was used

    if plotFig
        figure;
        hold on;
        xlim([0,PlotTimeRange]);
        
        legend('AutoUpdate','off');
        errorbar(time(1:40:PlotFrames),DistSqAvg(1:40:PlotFrames),Stdev(1:40:PlotFrames),'-.','CapSize',10,'LineWidth',1.2);
        AvgLegend = plot(time(1:PlotFrames),DistSqAvg(1:PlotFrames),'-','color','red','LineWidth',2);
        yrange = ylim;
        Height = yrange(2) - yrange(1);
        
        TimeFit = linspace(time(1),time(end),500);
        DistFit = polyval(coeff, TimeFit);
        plot(TimeFit,DistFit,'--','Color','black','LineWidth',1.5);
        
        for it = 1:N
            plot(Time{it}(1:PlotFrames),DistSq{it}(1:PlotFrames),'Color',[.7,.7,.7,.8])
        end

        if ShowText
            linFitEq = ['$r^{2} = ', sprintf('%.3g',coeff(1)),' t$'];
            RsqText = ['$ R^{2} = ',sprintf('%.4g',Rsq),'$'];
            text(40, yrange(1)+0.7*Height, linFitEq,'Interpreter','latex','FontSize',15);
            text(40, yrange(1)+0.8*Height, RsqText,'Interpreter','latex','FontSize',15);
            text(40, yrange(1)+0.9*Height,['$\eta = ',sprintf('%.5g',Viscosity),'$'],'Interpreter','latex','FontSize',15);
        end
              
        set(gca,'XScale','linear','YScale','linear','fontsize',20);
        legend(AvgLegend,'$\langle r^{2} \rangle$','Interpreter','latex','Location','northeast','FontSize',25);
        xlabel('$\mathrm{time} \left( s \right)$','Interpreter','latex','FontSize',25);
        ylabel('$\mathrm{distance^{2}} \left( \mu \mathrm{m}^{2} \right)$','Interpreter','latex','FontSize',25);
        title(['$\mathrm{Brownian \ motion} \left(', sprintf('%.15g',D),'\mu m \right)$'],'Interpreter','latex','FontSize',30);
        hold off;
    end
    
    switch nargout
        case 4
            varargout = {Time, X, Y, Viscosity};
        case 6
            varargout = {Time, X, Y, Viscosity, DistSqAvg, DistSqVar};
        case 7
            if CorrectDrift_default
                varargout = {Time, X, Y, Viscosity, DistSqAvg, DistSqVar, DriftVel};
            else
                error('ERR: ''DriftVel'' output can only be output when ''CorrectDrift'' option is used without input drift velocity');
            end
        otherwise
            error('ERR: number of output variables must be either 4, 6, or 7');
    end
    
end