clear; clc;  

% COHORTS
cohort(1).files     = pwd+"\Results\paper1\Raw data\ICU\MAT\";      % path to database
cohort(1).save      = pwd+"\Results\paper1\Processed data\ICU\";    % path to save folder
cohort(1).patients  = 140;                                          % number of patients in the cohort
cohort(1).duration  = 1200;                                         % length of analysis (in seconds) - "all" for whole recording
cohort(1).diff_pat  = [52,60,62,67,68,70,131,132,133,138];          % patients with different recording setting (reference and sampling frequency)
cohort(1).start     = 0;                                            % analysis start time (in seconds) - 0 for recording beginning

cohort(2).files     = pwd+"\Results\paper1\Raw data\SeLECTS\MAT\";
cohort(2).save      = pwd+"\Results\paper1\\Processed data\SeLECTS\";
cohort(2).patients  = 32;
cohort(2).duration  = 240;
cohort(2).diff_pat  = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28];
cohort(2).start     = 0;

% PARAMETERS
fbands      = [1,4;4,8;8,12;12,20;1,12;1,20];                       % frequency bands of interest
des_fs      = 250;                                                  % desired sampling frequency
channs      = ["F3","P3","Fz","Cz","F4","C4","P4","Fp1","Fp2","F7","T3","T5","O1","O2","T6","T4","F8","Pz"]; % EEG channels of interest
des_ref     = "C3";                                                 % desired reference
reref       = ["","average","mastoids","laplacian"];                % re-referencing techniques of interest [none,average,mastoids,laplacian]
rej_t       = 200;                                                  % artifact rejection threshold (µV)
window      = 2;                                                    % window  size (in seconds)
overlap     = window/2;                                             % overlap between windows (in seconds)
measrs      = ["CC","CORR","COH","IM","PLV","wPLI","MI","wSMI"];    % association measures of interest [CC,CORR,COH,IM,PLV,PLI,wPLI,MI,wSMI]
featrs      = 5;                                                    % number of connectivity features to extract

m           = measures(measrs);                                     % initialise measure object
m.mi.bins   = 300;                                                  % TO SPECIFY IF USING MI (put to 0 otherwise) - bin number 
m.wsmi.tau  = [80,40,28,16,28,16]./1000;                            % TO SPECIFY IF USING wSMI (put to 0 otherwise) - step tau
net         = network(featrs);                                      % initialise network object

% CONNECTIVITY
for c = 1:size(cohort,2)
    conn_mat = zeros(cohort(c).patients,length(reref),length(measrs),size(fbands,1),length(channs),length(channs));
    feat_mat = zeros([size(conn_mat,1,2,3,4),featrs]);

    for p = 1:cohort(c).patients
        disp("Processing patient n°" + p)
        pat_file = load(cohort(c).files + p + ".mat");
        eeg_data = data(pat_file.EEG,cohort(c).start,cohort(c).duration,des_fs,channs,des_ref); % initialise data object
        if ismember(p,cohort(c).diff_pat); eeg_data = eeg_data.original_ref; eeg_data = eeg_data.resampling(des_fs); end

        for r = 1:length(reref)
            % Pre-processing
            eeg_data = eeg_data.preprocess(reref(r),fbands,window,overlap,rej_t);
            eeg_data = eeg_data.add_ppsteps(m.pp_steps,fbands);
            
            % Associations computation
            ppdata_bands = eeg_data.ppdata_bands; ppdata_hann = eeg_data.ppdata_hann; fs = eeg_data.fs;
            for i = 1:length(channs)
                x_bands = []; x_hann = [];
                if ~isempty(ppdata_bands); x_bands = squeeze(ppdata_bands(:,i,:,:)); end
                if ~isempty(ppdata_hann); x_hann = squeeze(ppdata_hann(:,i,:)); end
                assoc = zeros(length(measrs),size(fbands,1),length(channs));
                parfor j = i+1:length(channs)
                    y_bands = []; y_hann = [];
                    if ~isempty(ppdata_bands); y_bands = squeeze(ppdata_bands(:,j,:,:)); end
                    if ~isempty(ppdata_hann); y_hann = squeeze(ppdata_hann(:,j,:)); end
                    assoc(:,:,j) = m.associations(x_bands,y_bands,x_hann,y_hann,fs,fbands);
                end
                conn_mat(p,r,:,:,i,:) = assoc(m.indicators_idx,:,:);
            end

            % Features extraction
            feat_mat(p,r,:,:,:) = net.compute_features(squeeze(conn_mat(p,r,:,:,:,:)),m);
        end
    end
    save(cohort(c).save + "conn_mat.mat", "conn_mat");
    save(cohort(c).save + "feat_mat.mat", "feat_mat");
end
