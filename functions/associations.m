classdef associations
    properties
        data            % data object for connectivity computation
        measure         % measure object for connectivity computation
        f_bands         % frequency bands of interest
        reref           % re-referencing schemes
        indicators      % association measures of interest
        indicators_idx  % indices of measures in the connectivity matrix
        window          % window size for connectivity computation
        pp_steps        % pre-processing steps to apply to each measure
        assoc_steps     % association steps to apply to each measure
        select_info     % choice of information measures
    end

    methods
        function obj = associations(eeg_data,measure,f_bands,reref,selected_measures,window)
            obj.data = eeg_data; obj.measure = measure;
            obj.f_bands = f_bands; obj.reref = reref; obj.window = window;
            [obj.indicators, obj.indicators_idx] = measure.get_measures(selected_measures);
            for i = 1:length(obj.indicators)
                obj.pp_steps = [obj.pp_steps, obj.indicators{i}.pp_steps];
                obj.assoc_steps = [obj.assoc_steps, obj.indicators{i}.assoc_steps];
            end
            obj.pp_steps = unique(obj.pp_steps,'stable');
            obj.assoc_steps = unique(obj.assoc_steps,'stable');
            if any(strcmp(selected_measures,"MI")); obj.select_info(1) = 1; end
            if any(strcmp(selected_measures,"wSMI")); obj.select_info(2) = 1; end
        end

        % Associations computation
        function conn_mat = compute_associations(obj)
            conn_mat = zeros(length(obj.reref),length(obj.indicators),size(obj.f_bands,1),length(obj.data.channels),length(obj.data.channels));
            for r = 1:length(obj.reref)
                % Pre-processing
                pp_data = obj.data.preprocess(obj.reref(r),obj.f_bands,obj.window);
                pp_data = pp_data.add_ppsteps(obj.pp_steps,obj.f_bands);
                % Association measures
                conn_mat(r,:,:,:,:) = obj.association_measures(pp_data);
            end
        end

        % Measures computation
        function conn_mat = association_measures(obj,pp_data)
            conn_mat = zeros(length(obj.indicators),size(obj.f_bands,1),length(obj.data.channels),length(obj.data.channels));
            parfor i = 1:length(obj.data.channels)
                x_bands = []; y_bands = []; x_hann = []; y_hann = [];
                if ~isempty(pp_data.ppdata_bands); x_bands = squeeze(pp_data.ppdata_bands(:,:,i,:)); end
                if ~isempty(pp_data.ppdata_hann); x_hann = squeeze(pp_data.ppdata_hann(:,i,:)); end
                assoc = zeros(length(obj.indicators),size(obj.f_bands,1),length(obj.data.channels));
                for j = i+1:length(obj.data.channels)
                    if ~isempty(pp_data.ppdata_bands); y_bands = squeeze(pp_data.ppdata_bands(:,:,j,:)); end
                    if ~isempty(pp_data.ppdata_hann); y_hann = squeeze(pp_data.ppdata_hann(:,j,:)); end
                    assoc(:,:,j) = obj.measure.assoc_steps(obj.assoc_steps,obj.indicators,x_bands,y_bands,x_hann,y_hann,pp_data.fs,obj.f_bands,obj.select_info);
                end
                conn_mat(:,:,i,:) = assoc(obj.indicators_idx,:,:,:);
            end
        end
    end
end
