classdef data
    properties
        samples     {mustBeNumeric} % number of data points for analysis
        channels                    % EEG channel names
        fs          {mustBeNumeric} % sampling frequency
        raw_data                    % EEG data [samples x channels]
        ref_data                    % EEG ref data [samples x 1]
        mastoids                    % EEG mastoids data [samples x 2]
        ppdata                      % pre-processed data (basic) [window size x channels x nb windows]
        ppdata_bands = [];          % pre-processed data (filtered in narrow bands) [bands x window size x channels x nb windows]
        ppdata_hann = [];           % pre-processed data (tapered by hanning window) [window size x channels x nb windows]
        padding_length = 0.5;       % padding size (in seconds) before narrow band filtering or hanning
        laplace = struct(...        % laplace re-referencing scheme
            'P3', ["T5","O1","Pz"], ...
            'T3', ["F7","T5","A1"], ...
            'T5', ["T3","P3","O1"], ...
            'F3', ["Fp1","F7","Fz"], ...
            'F7', ["Fp1","F3","T3"], ...
            'C3', ["F3","Cz","P3","T3"], ...
            'C4', ["F4","Cz","P4","T4"], ...
            'P4', ["C4","Pz","O2","T6"], ...
            'O1', ["T5","P3","Pz","O2"], ...
            'O2', ["T6","P4","Pz","O1"], ...
            'T6', ["T4","P4","O2","C4"], ...
            'T4', ["F8","C4","T6","A2"], ...
            'F4', ["Fp2","Fz","C4","F8"], ...
            'F8', ["Fp2","F4","T4","C4"], ...
            'Fp1', ["F7","F3","Fz","Fp2"], ...
            'Fp2', ["Fz","F4","F8","Fp1"], ...
            'Pz', ["O1","O2","P3","Cz","P4"], ...
            'Fz', ["Fp1","Fp2","F3","Cz","F4"], ...
            'Cz', ["Fz","Pz","C4","F3","F4","P4","P3"]);
    end

    methods
        function obj = data(EEG,start,duration,channels,ref)
            obj.fs = EEG.channels.Fs;
            if strcmp(duration,'all'); obj.samples = size(EEG.data,1);
            else; obj.samples = duration*obj.fs; end

            for ch = 1:size(EEG.data,2)
                chan = string(EEG.channels(ch).Name);
                if any(strcmp(chan, channels)) && ~strcmp(chan, ref)
                    obj.channels = horzcat(obj.channels, chan);
                    obj.raw_data = horzcat(obj.raw_data, EEG.data(start*obj.fs+1:end,ch));
                end
                if any(strcmp(chan, ["A1","A2"])); obj.mastoids = horzcat(obj.mastoids, EEG.data(start*obj.fs+1:end,ch)); end
                if any(strcmp(chan, ref)); obj.ref_data = EEG.data(start*obj.fs+1:end,ch); end
            end
        end

        % Adjustement of patients with different reference and/or sampling frequency
        function data = patient_to_homogenize(obj,des_fs)
            data = obj;
            % Re-referencing to desired reference
            data.raw_data = obj.raw_data - obj.ref_data;
            data.mastoids = obj.mastoids - obj.ref_data;
            % Downsampling to desired sampling frequency
            [p, q] = rat(des_fs/obj.fs);
            data.raw_data = resample(data.raw_data,p,q);
            data.mastoids = resample(data.mastoids,p,q);
            data.fs = des_fs;
        end

        % Data pre-processing
        function obj = preprocess(obj,reref,filter,window)
            reref_data = obj.rereferencing(reref);                                  % Re-referencing
            filt_data = obj.filter(reref_data,min(filter(:,1)),max(filter(:,2)));   % Broad band filtering
            obj.ppdata = obj.artefact_rejection(filt_data,window);                  % Artefacts rejection and windowing
        end

        % Additional pre-processing steps (association dependent)
        function obj = add_ppsteps(obj,pp_steps,f_bands)
            for step = pp_steps
                switch step
                    case "bands"; obj = obj.narrow_filtering(f_bands);              % Narrow band filtering
                    case "hann"; obj = obj.hanning_windows;                         % Hanning
                end
            end
        end

        % Re-referencing scheme choice
        function reref_data = rereferencing(obj,mode)
            if      strcmp(mode,"average");     reref_data = obj.rereference_to_average;
            elseif  strcmp(mode,"mastoids");    reref_data = obj.rereference_to_mastoids;
            elseif  strcmp(mode,"laplacian");   reref_data = obj.rereference_to_laplacian;
            else;                               reref_data = obj.raw_data; 
            end
        end

        % CAR re-referencing
        function reref_data = rereference_to_average(obj)
            ref = obj.raw_data;
            if any(strcmp(obj.channels,"Fp1")); ref(:,obj.channels=="Fp1")=[];  end
            if any(strcmp(obj.channels,"Fp2")); ref(:,obj.channels=="Fp2")=[];  end
            if any(strcmp(obj.channels,"Cz"));  ref(:,obj.channels=="Cz")=[];   end
            reref_data = obj.raw_data - mean(ref,2);
        end

        % Mastoids re-referencing
        function reref_data = rereference_to_mastoids(obj)
            reref_data = obj.raw_data - mean(obj.mastoids,2);
        end

        % Laplacian re-referencing
        function reref_data = rereference_to_laplacian(obj)
            reref_data = zeros(size(obj.raw_data));
            full_channels = [obj.channels,"A1","A2"];
            full_data = horzcat(obj.raw_data,obj.mastoids);
            for i = obj.channels
                if any(strcmp(i,["Fz","Pz"]))
                    weighted = mean(full_data(:,ismember(full_channels,obj.laplace.(i)(1:2))),2);
                    rest = sum(full_data(:,ismember(full_channels,obj.laplace.(i)(3:end))),2);
                    new_ref = (1/(length(obj.laplace.(i))-1))*(weighted+rest);
                else; new_ref = mean(full_data(:,ismember(full_channels,obj.laplace.(i))),2);
                end
                reref_data(:,obj.channels==i) = obj.raw_data(:,obj.channels==i) - new_ref;
            end
        end

        % Filtering
        function filtered_data = filter(obj,pp_data,f_highpass,f_lowpass)
            [b_high, a_high] = butter(1,f_highpass/(obj.fs/2),"high"); 
            [b_low, a_low] = butter(2,f_lowpass/(obj.fs/2));
            high_flit = filtfilt(b_high,a_high,pp_data);
            filtered_data = filtfilt(b_low,a_low,high_flit);
        end

        % Window-based artefacts rejection
        function wind_data = artefact_rejection(obj,pp_data,window)
            window = floor(window*obj.fs); overlap = floor(window/2);
            nb_wind = ceil((obj.samples-overlap)/(window-overlap));
            wind_data = []; w = 0;
            while size(wind_data,3) < nb_wind
                w_start = (window-overlap)*w+1;
                if w_start+window-1 > size(pp_data,1); warning("Too many artefacts. Analysis of "+size(wind_data,3)+"/"+nb_wind+" windows."); return; end
                wind = pp_data(w_start:w_start+window-1,:);
                if any(any(wind>200)) || any(any(wind<-200)) || any(max(wind)-min(wind)<2); w = w + 1;
                else; wind_data = cat(3,wind_data,wind); w = w + 1; 
                end
            end
        end

        % Symmetric window padding
        function [padded_data, padding_size] = padding_windows(obj,data)
            padding_size = obj.padding_length*obj.fs;
            pad_start = flip(data(1:padding_size,:,:));
            pad_end = flip(data(end-padding_size+1:end,:,:));
            padded_data = vertcat(pad_start,data,pad_end);
        end

        % Narrow bands filtering
        function obj = narrow_filtering(obj,f_bands)
            [extended_data, padding] = obj.padding_windows(obj.ppdata);
            filtered_data = zeros([size(f_bands,1),size(extended_data)]);
            for band = 1:size(f_bands,1); filtered_data(band,:,:,:) = obj.filter(extended_data,f_bands(band,1),f_bands(band,2)); end
            obj.ppdata_bands = filtered_data(:,padding+1:end-padding,:,:); 
        end

        % Normalization & Hanning tapering
        function obj = hanning_windows(obj)
            norm_data = (obj.ppdata - mean(obj.ppdata,1))./std(obj.ppdata,0,1);
            extended_data = obj.padding_windows(norm_data);
            hann_wind = repmat(hann(size(extended_data,1)),[1 size(extended_data,2,3)]);
            obj.ppdata_hann = hann_wind.*extended_data;
        end
    end
end
