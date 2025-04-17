clear;
lambda = 600;
I_0 = 2e4;
thetaShift = 0;
Ishift = 0;
Params = [lambda, I_0, thetaShift, Ishift];
slitDist = 200;
slitWidth = 85;
Beta = @(lambda, slitWidth, theta) 1e3*(pi/lambda)*slitWidth*sin(1e-3*theta);
Gamma = @(lambda, slitDist, theta) 1e3*(pi/lambda)*slitDist*sin(1e-3*theta);
theta = linspace(-10,10,100);
I = zeros(numel(theta),1);
for it = 1:numel(I)	
	if theta(it) ~= Params(3)
		I(it) = Params(2) * cos( Gamma(Params(1), slitDist, theta(it)-Params(3)) )^2;
		I(it) = I(it) * sin( Beta(Params(1), slitWidth, theta(it)-Params(3)) )^2 / Beta(Params(1), slitWidth, theta(it)-Params(3))^2;
		I(it) = I(it) + Params(4);
	else
		I(it) = Params(2) * cos( Gamma(Params(1), Params(3), theta(it)) )^2 + Params(4);% + 200*(rand-0.5);
	end
end
[gamma, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitIntensity',lambda,thetaShift,Ishift,'Gauss','-v')
%[gamma, lambda, I_0] = ConvFit(theta,I,slitDist,slitWidth,'FitLambda',thetaShift,Ishift,'Gauss','-v')