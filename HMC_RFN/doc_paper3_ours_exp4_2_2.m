%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Paper 3 (doc_paper3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Proposed Method (doc_paper3_ours)
% Multi-emitter interleaved scenario (exp4) -
% PRI modulation types of different sources are inconsistent - Special Scenario (exp4_2) - 
% Process in parallel between different Monte Carlo experiments (exp4_2_2)
close all; clear;

%-------------------------------------------------------------------------0
%------------------- Experiment Parameter Settings -------------
EmitterNum = 2;                                       % Number of emitters
Vnum = 5;                                             % Number of environmental variable changes
IterNum = 2000;                                       % Number of Monte Carlo trials under current conditions
DeleavingResult = [];
Time = zeros(Vnum, IterNum);                          % Storage for sorting time
T = 1e4;                                              % Total observation duration
Pm0 = 0.1;                                            % Pulse loss rate - target source parameter
Sp = 0.1;                                             % Spurious pulse ratio - target source parameter
PRIJitter = 0.005;                                    % PRI jitter rate - target source parameter
Sp_setting = 0.05;                                    % Spurious pulse ratio - RFN parameter
Pm_setting = 0.50;                                    % Pulse loss rate - RFN parameter
Alpha = 10;                                           % Maximum number of consecutive lost pulses - RFN parameter
Scenario_label = 2;                                   % PRI modulation parameters -- target pulse sequence parameters
%---------------------Performance Metrics------------------------------
Psearch1 = zeros(Vnum, IterNum);                       % Performance metric
Psearch2 = zeros(Vnum, IterNum);                       % Performance metric
Psearch3 = zeros(Vnum, IterNum);                       % Performance metric
%-------------------------------------------------------------------------1        
% Start MATLAB parallel pool
if isempty(gcp('nocreate'))
    parpool('local', 40);  % Adjust according to actual needs, e.g., parpool('local', N) where N is the number of threads
end

%----------------------------- Monte Carlo Experiment Loop ---------------------------
for sn = 1:Vnum 
    disp(['Iteration ' num2str(sn) ' (out of ' num2str(Vnum) ')']);
    Pm = Pm0 * (sn - 1);
    parfor j = 1:IterNum 
        disp(['Iteration ' num2str(sn) '/' num2str(Vnum) ', Step ' num2str(j) '/' num2str(IterNum)]);
        %-----------------------------------------------------------------0
        %---------------------- Interleaved TOA - Parameter Settings ----------------------------
        Input1 = cell(1, EmitterNum);
        Input2 = cell(1, EmitterNum);
        p_n_set = cell(1, EmitterNum);
        
        if Scenario_label == 1
            %------------------ Emitter 1
            Jagging_Order = randi([2, 5]);
            unique_floats = randperm(51, Jagging_Order) + 49 + (rand(1, Jagging_Order) - 0.5);
            Input1(1) = {'JaggingPri'};                                    % PRI modulation type -- Staggered
            Input2(1) = {unique_floats};                                   % PRI modulation parameters
            p_n_set{1} = unique_floats;
            %------------------ Emitter 2
            Input1(2) = {'FixedPri'};                                      % PRI modulation type -- FixedPri
            Input2(2) = {sum(unique_floats)};                              % PRI modulation parameters
            p_n_set{2} = sum(unique_floats);
            
        elseif Scenario_label == 2
            %------------------ Emitter 1
            DS_Order = randi([2, 5]);
            unique_floats = randperm(51, DS_Order) + 99 + (rand(1, DS_Order) - 0.5);
            unique_integers = randperm(7, DS_Order) + 1;
            Input1(1) = {'GroupPri'};                                      % PRI modulation type -- D&S
            Input2(1) = {[unique_floats; unique_integers]};                % PRI modulation parameters
            temp = [unique_floats; unique_integers];
            p_n = [];
            for n = 1:length(temp(1,:))
                p_n = [p_n, repmat(temp(1,n), 1, temp(2,n))];
            end
            p_n_set{1} = p_n;
            %------------------ Emitter 2
            Input1(2) = {'FixedPri'};                                      % PRI modulation type -- FixedPri
            Input2(2) = {unique_floats(1)};                                % PRI modulation parameters
            p_n_set{2} = unique_floats(1);
        end
        Input3 = repmat(1000, 1, EmitterNum);                              % Number of pulses (actual number of pulses constrained by the observation market)
        Input4 = repmat(0.001, 1, EmitterNum);                             % Time measurement error
        Input5 = repmat(Pm, 1, EmitterNum);                                % Pulse loss rate
        Input6 = 200 * rand(1, EmitterNum);                                % Start time
        Input7 = repmat(T, 1, EmitterNum);                                 % Observation duration
        Input8 = repmat(PRIJitter, 1, EmitterNum);                         % PRI jitter rate
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
        TOA_idx = 1:TOALen;
        %-----------------------------------------------------------------1
        
        %-----------------------------------------------------------------0
        %----------------------------- Pulse Sorting ------------------------------
        tic;                                          % Start timing
        Target = cell(1, EmitterNum);                                       % Preallocate Target array to store results
        for i = 1:EmitterNum
            % Construct Residual Fence Network
            [WO, OO, OE, WE] = ConstructionOfRFN(TOA, p_n_set{i}, Pm_setting, PRIJitter, Sp_setting, Alpha);
            % Function Purpose: Construct Residual Fence Network
            % Detailed description of inputs and outputs:
            % Inputs:
            % Input1: Interleaved pulse stream / Input2: PRI modulation information (target emitter source parameters) /
            % Input3: Set pulse loss rate (RFN hyperparameter) / Input4: PRI jitter rate (target emitter source parameters) /
            % Input5: Spurious pulse ratio (RFN hyperparameter) / Input6: Maximum consecutive lost pulse count (RFN hyperparameter)
            % Outputs:
            % Output1: Weights from starting point to intermediate nodes / Output2: Weights from intermediate nodes to intermediate nodes
            % Output3: Weights from intermediate nodes to endpoint / Output4: Weights from starting point to endpoint
            
            % Finding the optimal path based on the Residual Fence Network
            TOA_idx_temp = FindOptimalPathWithinRFN(WO, OO, OE, WE);
            % Function Purpose: Find the optimal path based on the Residual Fence Network
            % Detailed description of inputs and outputs:
            % Inputs:
            % Input1: Weights from starting point to intermediate nodes /
            % Input2: Weights from intermediate nodes to intermediate nodes /
            % Input3: Weights from intermediate nodes to endpoint /
            % Input4: Weights from starting point to endpoint /
            % outputArg1: Index of the pulses belonging to the target source
            
            TOA(TOA_idx_temp) = [];
            Target{i} = TOA_idx(TOA_idx_temp);
            TOA_idx(TOA_idx_temp) = [];
        end
        toc;                                          % End timing
        %-----------------------------------------------------------------1
        
        %-----------------------------------------------------------------0
        %----------------------------- Metric Calculation ------------------------------
        Psearch1_temp = zeros(1, EmitterNum);
        Psearch2_temp = zeros(1, EmitterNum);
        Psearch3_temp = zeros(1, EmitterNum);
        for i = 1:EmitterNum
            if ~isempty(Target{i})
                Index2 = zeros(1, TOALen);
                Index2(Target{i}) = 1;
                Correlation = sum(Index2 .* Index1(i, :));
                if Index1(i, end) == 0
                    Psearch1_temp(i) = Correlation.^2 / length(Target{i}) / cell2mat(PRILabel(i, 4));
                    Psearch2_temp(i) = Correlation / length(Target{i});
                    Psearch3_temp(i) = Correlation / cell2mat(PRILabel(i, 4));
                elseif Index1(i, end) == 1
                    Psearch1_temp(i) = Correlation.^2 / length(Target{i}) / (cell2mat(PRILabel(i, 4)) - 1);
                    Psearch2_temp(i) = Correlation / length(Target{i});
                    Psearch3_temp(i) = Correlation / (cell2mat(PRILabel(i, 4)) - 1);
                end
            else
                Psearch1_temp(i) = 0;
                Psearch2_temp(i) = 0;
                Psearch3_temp(i) = 0;     
            end
        end
        Psearch1(sn, j) = mean(Psearch1_temp);
        Psearch2(sn, j) = mean(Psearch2_temp);
        Psearch3(sn, j) = mean(Psearch3_temp);
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
%save('Mat.mat','Mat');
save(['doc_paper3_ours_exp4_2_2_', num2str(Scenario_label), '-', num2str(Sp_setting), '-', num2str(Pm_setting), '-', num2str(Sp), '-', num2str(Pm0), '-', num2str(Alpha), '-', num2str(EmitterNum), '.mat'], 'Mat');
