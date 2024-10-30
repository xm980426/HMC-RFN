function [Output1,Output2,Output3] = GenerateToaForOneWholeExp3(Input1,Input2,Input3,Input4,Input5,Input6,Input7,Input8,Input9,Input10,Input11)
% Generate interleaved pulse trains, including spurious pulses and considering jitter
% Inputs
% Input1: PRI type
% Input2: PRI parameters
% Input3: Number of pulses
% Input4: Time measurement error/us
% Input5: Pulse loss rate (<1)
% Input6: Start time/us
% Input7: Observation end time/us
% Input8: Jitter rate
% Input9: Spurious pulse ratio
% Input10: Current experimental scenario number
% Input11: Current Monte Carlo experiment iteration
% Outputs
% Output1: Interleaved TOA sequence
% Output2:  "Source ID, PRI modulation type, PRI modulation parameters,Number of pulse "
% Output3: Index of pulse from each source in the interleaved pulse sequence

EmitterNum = length(Input1);                                               % Number of sources  
PRILabel = cell(EmitterNum,4);                                             % Store source ID, PRI type, and its parameters
Toa = cell(1,EmitterNum);                                                  % Arrival time storage
Sp = Input9;                                                               % Spurious pulse ratio
VSn = Input10;                                                             % Current experimental scenario number
IterSn = Input11;                                                          % Current Monte Carlo experiment iteration
toaLen1 = zeros(1,EmitterNum);                                             % Store cumulative lengths
toaLen2 = zeros(1,EmitterNum);                                             % Store lengths of each segment
for k = 1:EmitterNum
    %% Generate PRI
    Toa(k) = {GenerateToaForPRIWithPRIJitter(Input1(k),Input2(k),Input3(k),Input4(k),Input5(k),Input6(k),Input7(k),Input8(k))};
    %% Generate TOA, GenerateToaForPRI is used to generate pulse arrival times for various repetition frequency patterns
    %% Inputs & Outputs
    % Inputs:
    % Input1: Repetition frequency type
    % Input2: Repetition frequency parameters
    % Input3: Number of pulses
    % Input4: Time measurement error/us
    % Input5: Pulse loss rate (<1)
    % Input6: Start time/us
    % Input7: Observation end time/us
    % Input8: Jitter rate
    % Outputs:
    % Output1: Pulse arrival times
    Len = length(cell2mat(Toa(k)));
    PRILabel(k,:) = [{[VSn,IterSn,k]}, Input1(k), Input2(k), Len];
    toaLen2(k) = Len;                                                      % Length of each segment
    toaLen1(k) = sum(toaLen2);                                             % Cumulative length
end
TOA = cell2mat(Toa);                                                       % Aggregate arrival times from all radiation sources
SpLen = round(length(TOA) * Sp / (1 - Sp));                                % Number of spurious pulses
if SpLen >= 2
    SpToa = max(TOA) * rand(1, SpLen - 2);                                 % Arrival times for spurious pulses
    TOA = [TOA, SpToa, max(TOA) + 100 * rand(1, 2)];                     
elseif SpLen == 1
    TOA = [TOA, max(TOA) + 100 * rand(1, 1)];                           
end
% Record the indices of each pulse sequence in the interleaved pulse train
[TOA, I] = sort(TOA);                                                       % Sort
Index = zeros(EmitterNum, length(TOA));                                     % Store indices of each radiation source pulse in the interleaved pulse train
Index(1,:) = (I <= toaLen1(1));
if EmitterNum >= 2
    for i = 2:EmitterNum
        Index(i,:) = (I > toaLen1(i-1)) .* (I <= toaLen1(i));             % Index of arrival times for each radiation source
    end
end
Output1 = TOA;
Output2 = PRILabel;
Output3 = Index;
end