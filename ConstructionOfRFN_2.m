function [outputArg1,outputArg2,outputArg3,outputArg4] = ConstructionOfRFN_2(inputArg1,inputArg2,inputArg3,inputArg4,inputArg5,inputArg6)
% Function Purpose: Construct Residual Fence Network
% Further subdivide the distribution range
% The PRI modulation type must be Fixed and the PRI jitter must be evenly distributed 
% Detailed description of inputs and outputs:
%% Inputs & Outputs
% Inputs:
% Input1: Interleaved pulse stream / Input2: PRI modulation information (target emitter source parameters) /
% Input3: Set pulse loss rate (RFN hyperparameter) / Input4: PRI jitter rate (target emitter source parameters) /
% Input5: Spurious pulse ratio (RFN hyperparameter) / Input6: Maximum consecutive lost pulse count (RFN hyperparameter)
% Outputs:
% Output1: Weights from starting point to intermediate nodes / Output2: Weights from intermediate nodes to intermediate nodes
% Output3: Weights from intermediate nodes to endpoint / Output4: Weights from starting point to endpoint 
%% ------------------------------------ Section Separator -------------------------------
% --------------------------------------------------------- Assign input parameters
TOA_interleaved = inputArg1;                                               % Interleaved pulse stream
p_n = inputArg2;                                                           % PRI modulation pattern
MissingRate = inputArg3;                                                   % Set pulse loss rate
JitterRate = inputArg4;                                                    % Set PRI jitter rate
SpuriousRate = inputArg5;                                                  % Set spurious pulse ratio
Alpha = inputArg6;                                                         % Maximum consecutive lost pulse count
%% ------------------------------------ Section Separator -------------------------------
% ---------------------------------------------------------- Construct Residual Fence Network
N = length(p_n);                                                           % Length of the state set
p_n = repmat(p_n,1,100);
M = length(TOA_interleaved)-1;                                             % Number of layers in the network
Y_n = diff(TOA_interleaved);                                               % Observation value vector
Belta_set = repmat((1-SpuriousRate)/N,1,N);                                % Initial feature vector
% -------------------------------------- 00: Calculate the state transition probability of the actual pulse sequence, f_ij
f_ij = zeros(N,N);                                                         % Declaration
for i = 1:N
    for j=1:N
        term = 0;
        if i<j
            u = 0;
            while(j-i-1+u*N<=Alpha)
                term = term + (1-MissingRate)*MissingRate.^(j-i-1+u*N);
                u = u+1;
            end
        elseif i>=j
            u = 1;
            while(j-i-1+u*N<=Alpha)
                term = term + (1-MissingRate)*MissingRate.^(j-i-1+u*N);
                u = u+1;
            end
        end
        f_ij(i,j) = term;
    end
end
% ------------------------------------ Section Separator --------------------------------
% ----------------------------------------- 02: Calculate the observation probability of the actual pulse sequence, d_ij
d_ij = zeros(N,Alpha+1);                                                   % Declaration
for i = 1:N
    for j = 1:Alpha+1
       d_ij(i,j) =  (1-MissingRate)*MissingRate.^(j-1);
    end
end
% ------------------------------------ Section Separator --------------------------------
% ------------------------------------------------- 03: Calculate weights from starting point to intermediate nodes 
Theta = M/(TOA_interleaved(end)-TOA_interleaved(1));                       % Intensity of Poisson flow
WO = zeros(M,N);                                                           % Declaration
for u = 1:M
    for i = 1:N
        if u ==1
            WO(u,i) = 1*Belta_set(i);
        elseif u >1
            y_set = Y_n(1:u-1);
            term = cumprod(exp(-1*y_set*Theta));
            WO(u,i) = term(end)*Belta_set(i)*SpuriousRate.^(u-1);
        end
    end
end
% ------------------------------------ Section Separator --------------------------------
% --------------------------------------------- 04: Calculate weights from intermediate nodes to intermediate nodes 
OO = zeros(M,N,M,N);
for u = 1:M
    for i = 1:N
        for v = 1:M
            for j = 1:N
                % Loop body -- Start
                if v<=u
                    OO(u,i,v,j) = 0;
                elseif v>u
                    % Secondary body -- Start
                    y_set = Y_n(u:v-1);
                    right = (1-SpuriousRate)*f_ij(i,j)*SpuriousRate.^(v-u-1);
                    if i<j
                        m = 0;                                             % Distinction
                        while(j-i-1+m*N<=Alpha)
                            p_n_set = p_n(i:j-1+m*N);
                            if sum(y_set)>=sum(p_n_set)*(1-JitterRate)&&sum(y_set)<=sum(p_n_set)*(1+JitterRate)
                                if v-u==1
                                    left = d_ij(i,j-i+m*N)/sum(d_ij(i,(j-i):N:(Alpha+1)));
                                elseif v-u>1
                                    term = cumprod( exp(-1*y_set(1:end-1)*Theta) );
                                    left = d_ij(i,j-i+m*N)/sum(d_ij(i,(j-i):N:(Alpha+1)))*term(end);
                                end
                                % The probability of the observation is adjusted according to the distribution probability of the subdivision in which the observation is located
                                OO(u,i,v,j) = Gain(sum(y_set),p_n_set,JitterRate)*left*right;
                                break;
                            else
                                m=m+1;
                            end
                        end
                    elseif i>=j
                        m = 1;                                             % Distinction
                        while(j-i-1+m*N<=Alpha)
                            p_n_set = p_n(i:j-1+m*N);
                            if sum(y_set)>=sum(p_n_set)*(1-JitterRate)&&sum(y_set)<=sum(p_n_set)*(1+JitterRate)
                                if v-u==1
                                    left = d_ij(i,j-i+m*N)/sum(d_ij(i,(j-i+N):N:(Alpha+1)));
                                elseif v-u>1
                                    term = cumprod(exp(-1*y_set(1:end-1)*Theta));
                                    left = d_ij(i,j-i+m*N)/sum(d_ij(i,(j-i+N):N:(Alpha+1)))*term(end);
                                end
                                % The probability of the observation is adjusted according to the distribution probability of the subdivision in which the observation is located
                                OO(u,i,v,j) = Gain(sum(y_set),p_n_set,JitterRate)*left*right;
                                break;
                            else
                                m=m+1;
                            end
                        end
                    end
                    % Secondary body -- End
                end
                % Loop body -- End
            end
        end
    end
end
% ------------------------------------ Section Separator --------------------------------
% ------------------------------------------------- 05: Calculate weights from intermediate nodes to endpoint 
OE = zeros(M,N);                                                           % Declaration
for u = 1:M
    for i = 1:N
        y_set = Y_n(u:M);
        term = cumprod(exp(-1*y_set*Theta));
        OE(u,i) = term(end)*SpuriousRate.^(M-u);
    end
end
% ------------------------------------ Section Separator --------------------------------
% ----------------------------------------------------- 06: Calculate weights from starting point to endpoint 
y_set = Y_n;
term = cumprod(exp(-1*y_set*Theta));
left = term(end);
right = SpuriousRate.^M;
WE = left*right;    
%% ------------------------------------ Section Separator -------------------------------
% ---------------------------------------------------- Output the results of the residual fence network construction 
outputArg1 = WO;
outputArg2 = OO;
outputArg3 = OE;
outputArg4 = WE;
end
