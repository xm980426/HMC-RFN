function [outputArg1] = FindOptimalPathWithinRFN(inputArg1,inputArg2,inputArg3,inputArg4)
% Function Purpose: Find the optimal path based on the Residual Fence Network
% Detailed description of inputs and outputs:
%% Inputs & Outputs
% Inputs:
% Input1: Weights from starting point to intermediate nodes /
% Input2: Weights from intermediate nodes to intermediate nodes /
% Input3: Weights from intermediate nodes to endpoint /
% Input4: Weights from starting point to endpoint /
% outputArg1: Index of the pulses belonging to the target source

%% ------------------------------------ Section Separator -------------------------------
% -------------------------------------------------------- Assign input parameters
WO = inputArg1;                                                            % Weights from starting point to intermediate nodes
OO = inputArg2;                                                            % Weights from intermediate nodes to intermediate nodes
OE = inputArg3;                                                            % Weights from intermediate nodes to endpoint
WE = inputArg4;                                                            % Weights from starting point to endpoint
%% ------------------------------------ Section Separator -------------------------------
% ----------------------------------------------- Find the optimal path based on the Residual Fence Network
SIZE = size(OO);
M = SIZE(1);
N = SIZE(2);
mul = 10;                                                                  % Very important to prevent weights from becoming zero after multiplication
% ------------------------------------ Section Separator --------------------------------
% ----------------------------------------------------------------- 00: Initialization 
phi = zeros(M,N);

psi = nan(M,N,2);

phi(1,:) = WO(1,:);
% ------------------------------------ Section Separator --------------------------------
% ------------------------------------------------------------------- 01: Induction 
for n = 1:M-1
    for i=1:N
        term = phi(1:n,:).*OO(1:n,:,n+1,i);
        [phi(n+1,i), linear_idx] = max(term(:));                           % Find the maximum element and its linear index
        
        if phi(n+1,i)>WO(n+1,i)
            [a,b] = ind2sub(size(term), linear_idx);                       % Convert linear index to row and column indices
            psi(n+1,i,1) = a; 
            psi(n+1,i,2) = b;   
        else
            phi(n+1,i) = WO(n+1,i);
            psi(n+1,i,1)=NaN;
            psi(n+1,i,2)=NaN;
        end
    end
    WO = WO.*mul;
    WE = WE*mul;
    phi(1:n+1,:) = phi(1:n+1,:).*mul;
end
% ------------------------------------ Section Separator --------------------------------
% ------------------------------------------------------------------- 02: Termination 
term = phi.*OE;
[phi_final, linear_idx] = max(term(:));                                    % Find the maximum element and its linear index
if phi_final>WE
    [a,b] = ind2sub(size(term), linear_idx);                               % Convert linear index to row and column indices
    psi_final(1) = a;
    psi_final(2) = b;
else
    psi_final=[NaN, NaN];
end
% --------------------------------------------------------------- 03: Path Backtracking 
x_n = zeros(1,M);
if N~=1
    while( ~isnan(psi_final(1)))
        x_n(psi_final(1)) = psi_final(2);
        psi_final(1:2) = psi(psi_final(1),psi_final(2),1:2);
    end
elseif N==1
     while( ~isnan(psi_final(1)))
        x_n(psi_final(1)) = 1;
        psi_final = psi(psi_final(1),1,1);
    end
end
%% ------------------------------------ Section Separator -------------------------------
% -------------------------------------------- Output the state sequence prediction results based on the optimal path 
outputArg1 = find(x_n>0);
end
