function metrics = computeDetectionMetrics(gt, alarms)
% This function computes detection performance metrics
% for the KEFSD Detector
% Inputs:
%       gt      : Nx1 ground truth vector (0 = normal, 1 = attack) where 
%                 N = no. of samples
%       alarms  : NxM matrix of alarm signals (0 or 1) where 
%                 M = number of detectors 
%
% Output:
%       TPR :      (1xM) True Positive Rate (Recall)
%       FAR :      (1xM) False Alarm Rate (False Positive Rate)
%       FNR :      (1xM) False Negative Rate
%       Precision :(1xM)
%       F1        :(1xM)

% Andrew Omiloli
% 03/01/2026
%% Ensure column vector
  gt = gt(:);

%% Check dimensions
if size(alarms,1) ~= length(gt)
  error('Ground truth and alarms must have the same number of samples');
end

%% Get number of detectors 
num_detectors = size(alarms,2);

%% Initialize outputs
TPR       = zeros(1, num_detectors);
FAR       = zeros(1, num_detectors);
FNR       = zeros(1, num_detectors);
Precision = zeros(1, num_detectors);
F1        = zeros(1, num_detectors);

%% Logical indices
attack_idx   = (gt == 1);
normal_idx   = (gt == 0);
total_attack = sum(attack_idx);
total_normal = sum(normal_idx);

%% Check content of ground truth
if total_attack == 0 || total_normal == 0
    error('Ground truth must contain both attack and normal samples');
end

for i = 1:num_detectors
    alarm = alarms(:,i);

    % Confusion matrix components
    TP = sum((alarm == 1) & attack_idx);
    FP = sum((alarm == 1) & normal_idx);
    % TN = sum((alarm == 0) & normal_idx);
    FN = sum((alarm == 0) & attack_idx);

    % Compute Metrics
    TPR(i) = TP / total_attack;        % Recall
    FAR(i) = FP / total_normal;        % False Positive Rate
    FNR(i) = FN / total_attack;

    %% When computing Precision, ensure TP + FP is not 0
    if (TP + FP) == 0
        Precision(i) = 0;
    else
        Precision(i) = TP / (TP + FP);
    end

    %% When computing F1, ensure Precision(i) + TPR(i) is not 0
    if (Precision(i) + TPR(i)) == 0
        F1(i) = 0;
    else
        F1(i) = 2 * (Precision(i) * TPR(i)) / (Precision(i) + TPR(i));
    end
end

    %% Store results
    metrics.TPR = TPR;
    metrics.FAR = FAR;
    metrics.FNR = FNR;
    metrics.Precision = Precision;
    metrics.F1 = F1;

end