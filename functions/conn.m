clear; clc;  

% PARAMETERS
f_bands     = [1,4;4,8;8,12;12,20;1,12;1,20];                       % frequency bands of interest
channels    = ["F3","P3","Fz","Cz","F4","C4","P4","Fp1","Fp2",...
                "F7","T3","T5","O1","O2","T6","T4","F8","Pz"];      % EEG channels of interest
ref         = "C3";                                                 % recording reference
fs          = 250;                                                  % sampling frequency
reref       = ["","average","mastoids","laplacian"];                % re-referencing techniques of interest
select_meas = ["CC","CORR","COH","IM","PLV","wPLI","MI","wSMI"];    % association measures of interest
window      = 2;                                                    % window  size (in seconds)
nb_features = 5;                                                    % number of connectivity features to extract

% COHORTS
cohort(1).files     = pwd + "\Results\paper\Raw data\ICU\MAT\";     % path to database
cohort(1).save      = pwd + "\Results\paper\Processed data\ICU\";   % path to save folder
cohort(1).start     = 0;                                            % analysis start time (in seconds) - 0 for recording beginning
cohort(1).duration  = 1200;                                         % length of analysis (in seconds) - "all" for whole recording
cohort(1).patients  = 140;                                          % number of patients in the cohort
cohort(1).diff_pat  = [52,60,62,67,68,70,131,132,133,138];          % patients with different recording setting

cohort(2).files     = pwd + "\Results\paper\Raw data\SeLECTS\MAT\";
cohort(2).save      = pwd + "\Results\paper\Processed data\SeLECTS\";
cohort(2).start     = 0;
cohort(2).duration  = 240;
cohort(2).patients  = 32;
cohort(2).diff_pat  = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28];

% CONNECTIVITY
measure = measures;
net = network(nb_features);
for c = 1:size(cohort,2)
    conn_mat = zeros(cohort(c).patients,length(reref),length(select_meas),size(f_bands,1),length(channels),length(channels));
    feat_mat = zeros([size(conn_mat,1,2,3,4),nb_features]);
    
    for p = 1:cohort(c).patients
        disp("Processing patient n°" + p)
        pat_file = load(cohort(c).files + p + ".mat");
        eeg_data = data(pat_file.EEG,cohort(c).start,cohort(c).duration,channels,ref);
        if ismember(p,cohort(c).diff_pat); eeg_data = eeg_data.patients_to_homogenize(fs); end

        % Connectivity matrices computations
        assoc = associations(eeg_data,measure,f_bands,reref,select_meas,window);
        matrix = assoc.compute_associations;
        conn_mat(p,:,:,:,:,:) = matrix;

        % Features extraction
        feat_mat(p,:,:,:,:) = net.compute_features(matrix,assoc);
    end
    save(cohort(c).save + "conn_mat.mat", "conn_mat");
    save(cohort(c).save + "feat_mat.mat", "feat_mat");
end
