function [outputArg1] = ProbabilityoftheSubdivisionInterval(inputArg1,inputArg2,inputArg3,inputArg4)
% Function Purpose: Function to calculate the distribution probability of the observation value in the subdivision interval
% The PRI modulation type must be Fixed and the PRI jitter must be evenly distributed 
% Detailed description of inputs and outputs:
%% Inputs & Outputs
% Inputs:
% Input1: Observation value / Input2: Number of uniform distribution variables in the joint distribution - 1
% Input3: Coefficient before uniform distribution variable (value of fixed PRI) / Input4: Distribution radius of uniform distribution variable (jitter rate)
% Outputs:
% Output1: Probability of the subdivision interval
x = inputArg1;
j = inputArg2;
p = inputArg3;
JitterRate = inputArg4;
if j == 0
    outputArg1 = 1;
else
    n = j+1;
    y = (x/p-n*(1-JitterRate))/(2*JitterRate);                             % Observation value of the joint random variable of (j+1) standard uniformly distributed random variables
    delta = 0.01;                                                          % Quantization interval                     
    idx = floor(y/0.5);                                                    % Index of the subdivision interval where the observation value is located
    z = idx*0.5:delta:(idx+1)*0.5-delta;                                   % The number of subdivision intervals is 2n
    k = 0:1:n;
    zz = repmat(z,length(k),1);
    kk = repmat(k.',1,length(z));
    % Calculate the distribution probability of the observation value in the subdivision interval through integration
    outputArg1 = 1/2./factorial(n-1)*sum(sum( (zz-kk).^(n-1) .* sign(zz-kk) .* (-1).^kk .* arrayfun(@(q) nchoosek(n, q), kk) ))*delta;
end
end
