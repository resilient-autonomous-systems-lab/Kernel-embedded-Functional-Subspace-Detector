%% ROC curve and AUC to measure detector performance 
% Note: Ensure detector output and attacks ground truth signals are 
% logged first! 
detector_signal = squeeze(out.logsout.get('detector_output').Values.Data); % detector signals per area
ground_truth    = squeeze(out.logsout.get('ground_truth').Values.Data);     % ground truth

figure
hold on 
for i = 1:size(detector_signal,2)
     scores_i = detector_signal(:,i);   % one detector
     ground_truth_i = ground_truth(:,i);

     [PD_i, PFA_i, ~, AUC_i] = perfcurve(ground_truth_i, scores_i, 1);

     plot(PFA_i, PD_i, 'DisplayName', ...
         sprintf('Det %d (AUC=%.3f)', i, AUC_i))
end
 xlabel('False Alarm Rate')
 ylabel('Detection Probability')
 title('ROC per Detector (PCA)')
legend show
grid on

%% Display performancee metrics
alarms = squeeze(out.logsout.get('alarm').Values.Data);
metrics = computeDetectionMetrics(ground_truth, alarms)