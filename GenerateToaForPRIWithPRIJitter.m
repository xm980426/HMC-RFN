function [Output1] = GenerateToaForPRIWithPRIJitter(Input1,Input2,Input3,Input4,Input5,Input6,Input7,Input8)
% This function generates pulse arrival times for various repetition frequency patterns
% The function's purpose is to statistically determine possible PRIs
%% Inputs & Outputs
% Inputs:
% Input1: Repetition frequency type
% Input2: Repetition frequency parameters
% Input3: Number of pulses
% Input4: Time measurement error/us
% Input5: Pulse loss rate (<1)
% Input6: Start time/us
% Input7: Observation time/us
% Input8: Jitter rate
% Outputs:
% Output1: Pulse arrival times

%% Input
Type = cell2mat(Input1);
Parameter = cell2mat(Input2);
Num = Input3;
Error = Input4;
LossRate = Input5;
StartTime = Input6;
EndTime = Input7;
PRIJitterRate = Input8;

%% Main Process
%---------------------------------------------------------------- Jitter distribution
% Truncated normal distribution
NormPd = makedist('Normal','mu',0,'sigma',PRIJitterRate/3);
TruncatePdForPRIJitter = truncate(NormPd,-PRIJitterRate,PRIJitterRate);
NormrndValueForPRIJitter = random(TruncatePdForPRIJitter,[1,Num]);

% Uniform distribution
%NormrndValueForPRIJitter = 2*(rand(1,Num)-0.5)*PRIJitterRate;
%-------------------------------------

toa = zeros(1,Num);                                                        % Define TOA sequence
toa(1) = StartTime;                                                        % Confirm the first pulse
if strcmp(Type,'FixedPri')                                                 % Fixed                                             
    % Output1 = GenerateToaForPRI('FixedPri',{[230]},500,0.5,0,9);figure;stem( Output1(1:20),ones(1,20));diff(Output1);
    pri = Parameter;
    for i = 1:Num-1
        toa(i+1) = pri*(1+NormrndValueForPRIJitter(i))+toa(i); 
    end
    
elseif strcmp(Type,'JaggingPri')                                           % Staggered
    % Output1 = GenerateToaForPRI('JaggingPri',{[230,290,370]},1000,0.5,0,9);figure;stem( Output1(1:20),ones(1,20));diff(Output1);
    pri = Parameter;
    for i = 1:Num-1
        toa(i+1) = pri(mod(i-1,length(pri))+1)*(1+NormrndValueForPRIJitter(i))+toa(i); 
    end
elseif strcmp(Type,'SlipPri')                                              % Sliding
    % Output1 = GenerateToaForPRI('SlipPri',[{[230,300]},{10}],1000,0.5,0,9);figure;plot( Output1,ones(1,length( Output1)));diff(Output1);
    pri1 = Parameter(1);
    pri2 = Parameter(2);
    NumofSlip = Parameter(3);
    pri = linspace(pri1,pri2,NumofSlip);
    for i = 1:Num-1
        toa(i+1) = pri(mod(i-1,length(pri))+1)*(1+NormrndValueForPRIJitter(i))+toa(i); 
    end
 
elseif strcmp(Type,'GroupPri')                                             % D&S
    % Output1 = GenerateToaForPRI('GroupPri',[{[230,300,410]},{[100,200,300]}],1000,0.5,0,9);figure;plot( Output1,ones(1,length( Output1))); diff(Output1);
    pri = Parameter(1,:);
    num = Parameter(2,:);
    PRI = [];
    for i = 1:length(pri)
        PRI = [PRI,repmat(pri(i),1,num(i))];
    end
    pri = PRI;
    for i = 1:Num-1
        toa(i+1) = pri(mod(i-1,length(pri))+1)*(1+NormrndValueForPRIJitter(i))+toa(i); 
    end
    
    
elseif strcmp(Type,'SinedPri')                                             % Wobulated
    % Output1 = GenerateToaForPRI('SinedPri',[{[230]},{[20]},{[0.0005]},{[0]}],1000,0.5,0,9);figure;plot( Output1,ones(1,length( Output1))); diff(Output1);
    Pri = Parameter(1);                                                    % Mean
    priAm = Parameter(2);                                                  % Amplitude
    Nperiod = Parameter(3);                                                % Modulation period
    Phi = Parameter(4);                                                    % Initial phase/rad
    pri = Pri+priAm*sin(2*pi/Nperiod.*(1:Num)+Phi);
    for i = 1:Num-1
        toa(i+1) = pri(mod(i-1,length(pri))+1)*(1+NormrndValueForPRIJitter(i))+toa(i); 
    end
end
toa = toa+unifrnd(-1,1,[1,length(toa)])*Error;                             % Add measurement error
limit=rand(1,length(toa));                                                 % Introduce lost
toa (limit<LossRate) = [];
%% Output
Output1 = toa(toa<=EndTime);
end