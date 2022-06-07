% plot tension vs time averaging across many runs.
% put this script and its dependencies in folders filled
% with prom_ mat files, which are the outputs of the ring
% simulation

tens_mean_store = [];
tens_sd_store = [];
tanalysisRange = 1:100;
for tanalysis = tanalysisRange
    tens_store = [];
	% edit the index range accordingly
	for ianalysis = 1:5;
        filename = ['prom_',num2str(ianalysis),'_',num2str(tanalysis),'min.mat'];
        if ~exist(filename,'file')
            continue
        else
            load(filename)
        end
        tension_circ;
        tens_store = [tens_store,mean(tens)];
    end
    tens_mean_store = [tens_mean_store, mean(tens_store)];
    tens_sd_store = [tens_sd_store, std(tens_store)];
end
errorbar(tanalysisRange,tens_mean_store,tens_sd_store)