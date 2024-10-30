%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Paper 3 (doc_paper3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Proposed Method (doc_paper3_ours)
% Incorporate the probabilistic distribution characteristics of observations into the consideration of observational probabilities - (exp5)
% Set PRI jitter rate and spurious pulse ratio of the target source as an environmental variable (exp5_2) - 
% Process in parallel between different Monte Carlo experiments (exp5_2_2)

close all; clear;

%-------------------------------------------------------------------------0
%------------------------------ Experiment Parameter Settings-------------------------------
Vnum2 = 6;                                                                 % Condition variable 2 - spurious pulse ratio
Vnum3 = 5;                                                                 % Condition variable 3 - PRI jitter rate
Vnum = Vnum2 * Vnum3;                                                      % Number of condition variables

Sp_set = (1:Vnum2) * 0.1 + 0.1;
JitterRate_set = ((1:Vnum3) - 1) * 0.05 + 0.001;
[A, B] = meshgrid(Sp_set, JitterRate_set);
coords = [A(:) B(:)];

IterNum = 1000;                                                            % Number of Monte Carlo trials under current conditions
DeleavingResult = [];
Time = zeros(Vnum, IterNum);                                               % Storage for sorting time
T = 1e4;                                                                   % Total observation duration
Pm = 0.4;                                                                  % Pulse loss rate - target source parameter
Sp_setting = 0.05;                                                         % Spurious pulse ratio - RFN parameter
Pm_setting = 0.50;                                                         % Pulse loss rate - RFN parameter
Alpha = 10;                                                                % Maximum number of consecutive lost pulses - RFN parameter
Yes_label = 0;                                                             % Whether to incorporate the probabilistic distribution characteristics of observations into the consideration of observational probabilities (1: yes; 0: no)
%------------------------------------------------------------------- Performance Metrics
Psearch1 = zeros(Vnum, IterNum);                                          % Performance metric
Psearch2 = zeros(Vnum, IterNum);                                          % Performance metric
Psearch3 = zeros(Vnum, IterNum);                                          % Performance metric

%-------------------------------------------------------------------------1
        
% Start MATLAB parallel pool
if isempty(gcp('nocreate'))
    parpool('local', 40);  % Adjust according to actual needs, e.g., parpool('local', N) where N is the number of threads
end

%----------------------------- Monte Carlo Experiment Loop ---------------------------
for sn = 1:Vnum 
    disp(['Iteration ' num2str(sn) ' (out of ' num2str(Vnum) ')']);
    % Environmental variable values
    Sp = coords(sn, 1);                                                     % Spurious pulse ratio - target source parameter
    JitterRate = coords(sn, 2);                                            % Jitter rate - target source parameter
    parfor j = 1:IterNum  
    %for j = 1:IterNum 
        disp(['Iteration ' num2str(sn) '/' num2str(Vnum) ', Step ' num2str(j) '/' num2str(IterNum)]);
        %-----------------------------------------------------------------0
        %---------------------- Interleaved TOA - Parameter Settings -----------------------------
        PRItype = randi(5,1,1);                                            % PRI modulation information - target source parameter
        if PRItype == 1
            %------------------ Fixed
            random_float = 100 + (160 - 100) * rand;
            Input1 = {'FixedPri'};                                         % PRI modulation type
            Input2 = {random_float};                                       % PRI modulation parameters
            p_n = cell2mat(Input2);                                        % PRI modulation information
        elseif PRItype == 2
            %------------------ Staggered
            Jagging_Order = randi([2, 6]);
            unique_floats = randperm(51, Jagging_Order) + 99 + (rand(1, Jagging_Order) - 0.5);           
            Input1 = {'JaggingPri'};                                       % PRI modulation type
            Input2 = {unique_floats};                                      % PRI modulation parameters
            p_n = cell2mat(Input2);
        elseif PRItype == 3
            %------------------ D&S
            DS_Order = randi([2, 6]);
            unique_floats = randperm(51, DS_Order) + 99 + (rand(1, DS_Order) - 0.5);      
            unique_integers = randperm(9, DS_Order) + 1;
            Input1 = {'GroupPri'};                                         % PRI modulation type
            Input2 = {[unique_floats; unique_integers]};                   % PRI modulation parameters
            temp = cell2mat(Input2);
            p_n = [];
            for n = 1:length(temp(1, :))
                p_n = [p_n, repmat(temp(1, n), 1, temp(2, n))];
            end
        elseif PRItype == 4
            %------------------ Sliding
            random_integers = randi([130, 160]);
            random_integers_2 = randi([5, 15]);
            Input1 = {'SlipPri'};                                          % PRI modulation type
            Input2 = {[100, random_integers, random_integers_2]};          % PRI modulation parameters
            temp = cell2mat(Input2);
            p_n = linspace(temp(1), temp(2), temp(3));
        elseif PRItype == 5
            %------------------ Wobulated
            random_float = 100 + (160 - 100) * rand;
            random_float_2 = 10 + (30 - 10) * rand;
            random_integers = randi([10, 20]);
            Input1 = {'SinedPri'};                                         % PRI modulation type
            Input2 = {[random_float, random_float_2, random_integers, 0]}; % PRI modulation parameters
            temp = cell2mat(Input2);
            p_n = temp(1) + temp(2) * sin(2 * pi / temp(3) .* (1:temp(3) + temp(4)));
         end
        Input3 = 1000;                                 % Number of pulses
        Input4 = 0.001;                                % Time measurement error
        Input5 = Pm;                                   % Pulse loss rate
        Input6 = 100 * rand(1, 1);                     % Start time
        Input7 = T;                                    % Observation end time
        Input8 = JitterRate;                           % Unintentional jitter rate
        %-----------------------------------------------------------------1
        
        %-----------------------------------------------------------------0
        %----------------------- Interleaved TOA - Generate TOA ----------------------------
        [TOA, PRILabel, Index1] = GenerateToaForOneWholeExp3(Input1, Input2, Input3, Input4, Input5, Input6, Input7, Input8, Sp, sn, j);
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
        % Output2: "Source ID, PRI modulation type, PRI modulation parameters, Number of pulses"
        % Output3: Index of the pulse from each source in the interleaved pulse sequence
        
        TOALen = length(TOA);                         % Total length of interleaved pulse trains
        %-----------------------------------------------------------------1
        
        %-----------------------------------------------------------------0
        %----------------------------- Pulse Sorting ------------------------------
        tic;                                                               % Start timing
        % Construct Residual Fence Network 
        if Yes_label == 0
            [WO, OO, OE, WE] = ConstructionOfRFN(TOA, p_n, Pm_setting, JitterRate, Sp_setting, Alpha);
        elseif Yes_label == 1
            [WO, OO, OE, WE] = ConstructionOfRFN_2(TOA, p_n, Pm_setting, JitterRate, Sp_setting, Alpha);
        end
        % Function Purpose: Construct Residual Fence Network
        % Detailed description of inputs and outputs:
        % Inputs:
        % Input1: Interleaved pulse stream / Input2: PRI modulation information (target emitter source parameters) /
        % Input3: Set pulse loss rate (RFN hyperparameter) / Input4: PRI jitter rate (target emitter source parameters) /
        % Input5: Spurious pulse ratio (RFN hyperparameter) / Input6: Maximum consecutive lost pulse count (RFN hyperparameter)
        % Outputs:
        % Output1: Weights from starting point to intermediate nodes / Output2: Weights from intermediate nodes to intermediate nodes
        % Output3: Weights from intermediate nodes to endpoint / Output4: Weights from starting point to endpoint
        
        % Find optimal path based on Residual Fence Network
        [TOASerialNumSingle] = FindOptimalPathWithinRFN(WO, OO, OE, WE);
        % Function Purpose: Find the optimal path based on the Residual Fence Network
        % Detailed description of inputs and outputs:
        % Inputs:
        % Input1: Weights from starting point to intermediate nodes /
        % Input2: Weights from intermediate nodes to intermediate nodes /
        % Input3: Weights from intermediate nodes to endpoint /
        % Input4: Weights from starting point to endpoint /
        % outputArg1: Index of the pulses belonging to the target source
        toc;                                                               % End timing
        %-----------------------------------------------------------------1
        
        %-----------------------------------------------------------------0
        %----------------------------- Metric Calculation ------------------------------
        TOASerialNumTemp = TOASerialNumSingle;
        if ~isempty(TOASerialNumTemp)
            Index2 = zeros(1, TOALen);
            Index2(TOASerialNumTemp) = 1;
            Correlation = sum(Index2 .* Index1(1, :));
             if Index1(1, end) == 0
                Psearch1(sn, j) = Correlation.^2 / length(TOASerialNumTemp) / cell2mat(PRILabel(1, 4));
                Psearch2(sn, j) = Correlation / length(TOASerialNumTemp);
                Psearch3(sn, j) = Correlation / cell2mat(PRILabel(1, 4));
            elseif Index1(1, end) == 1
                Psearch1(sn, j) = Correlation.^2 / length(TOASerialNumTemp) / (cell2mat(PRILabel(1, 4)) - 1);
                Psearch2(sn, j) = Correlation / length(TOASerialNumTemp);
                Psearch3(sn, j) = Correlation / (cell2mat(PRILabel(1, 4)) - 1);
             end
        else
            Psearch1(sn, j) = 0;
            Psearch2(sn, j) = 0;
            Psearch3(sn, j) = 0;
        end
        Time(sn, j) = toc;
        %-----------------------------------------------------------------1
    end

end


%% -------------------------- Save Results ------------------------------
Indic(:, 1) = mean(Psearch1, 2);
Indic(:, 2) = mean(Psearch2, 2);
Indic(:, 3) = mean(Psearch3, 2);
Indic(:, 4) = mean(Time, 2);
%%
Mat = Indic;
save(['doc_paper3_ours_exp5_2_2_', num2str(Yes_label), '-', num2str(Sp_setting), '-', num2str(Pm_setting), '-', num2str(Alpha), '-', num2str(Pm), '.mat'], 'Mat');
