function [outputArg1] = Gain(inputArg1, inputArg2, inputArg3)
% Function purpose: returns the gain determined by the probability density function
% Detailed description of inputs and outputs:

x = inputArg1;                                                             % Observation value
p_n_set = inputArg2;                                                       % p_n sequence
JitterRate = inputArg3;                                                    % Jitter rate
mu = sum(p_n_set);                                                         % Mean
sigma = sqrt(sum(p_n_set.^2)) * JitterRate / 3;                            % Standard deviation
Min = mu * (1 - JitterRate);                                               % Lower bound of the truncated interval
Max = mu * (1 + JitterRate);                                               % Upper bound of the truncated interval

% Determine the interval to which the observation value belongs
NumberOfInterval = 3 * length(p_n_set);                                    % Number of subdivided intervals
Edge = linspace(Min, Max, NumberOfInterval + 1);                            % Boundaries of the subdivided intervals
Idx = discretize(x, Edge);                                                % Subdivided interval to which the observation value belongs

% Calculate the distribution probability of the interval [a, b]
% probability = normcdf(b, mu, sigma) - normcdf(a, mu, sigma);
probability = normcdf(Edge(Idx + 1), mu, sigma) - normcdf(Edge(Idx), mu, sigma);
% Calculate gain based on the distribution probability
gain = probability / (1 / NumberOfInterval);
% Output result
outputArg1 = gain;
end
