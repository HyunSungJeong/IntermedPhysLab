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

    %{
    Params0 = [1e-5, 600, 1500];
    [gamma, lambda, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitLambda',4.4,0,'Gauss','InitParams',Params0,'-v');
    disp(['gamma = ',sprintf('%.15g',gamma)]);
    disp([lambda, I_0]);
    %}

    %{
    Params0 = [600, 1500];
    [lambda, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitLambda',4.4,0,'InitParams',Params0,'-v');
    disp([lambda, I_0]);
    %}

    %{}
    Params0 = [1e-5, 1500];
    [gamma, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitIntensity',670,4.4,0,'Lorentz','InitParams',Params0,'-v');
    disp(['gamma = ',sprintf('%.15g',gamma)]);
    disp(I_0);
    %}
end