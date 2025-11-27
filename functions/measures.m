classdef measures
    properties
        all_names = ["CC","CORR","COH","IM","PLV","PLI","wPLI","MI","wSMI"];
        cc      = struct('name','cc','pp_steps',"bands",'assoc_steps',"cc",'function',@computeCC,'max_lag',0.25);  
        corr    = struct('name','corr','pp_steps',"bands",'assoc_steps',"cc",'function',@computeCORR); 
        coh     = struct('name','coh','pp_steps',"hann", 'assoc_steps',["cpsd","cohy"],'function',@computeCOH); 
        im      = struct('name','im','pp_steps',"hann",'assoc_steps',["cpsd","cohy"],'function',@computeIM);
        plv     = struct('name','plv','pp_steps',"hann",'assoc_steps',"cpsd",'function',@computePLV);      
        pli     = struct('name','pli','pp_steps',"hann",'assoc_steps',"cpsd",'function',@computePLI);     
        wpli    = struct('name','wpli','pp_steps',"hann",'assoc_steps',"cpsd",'function',@computeWPLI);   
        mi      = struct('name','mi','pp_steps',"bands",'assoc_steps',"info",'function',@computeMI,'bins',300);                            
        wsmi    = struct('name','wsmi','pp_steps',"bands",'assoc_steps',"info",'function',@computeWSMI,...
                    'kernel',3,'tau',[80,40,28,16,28,16]./1000,'weights',[0,1,1,1,0,1;1,0,1,1,1,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,1,0,1;1,0,1,1,1,0]);

        f_bands     % Frequency bands of interest    
        cross_corr  % Cross-correlogram
        cross_spec  % Cross-spectrum
        cohy        % Coherency
        wind_mi     % Mutual information on windows
        wind_wsmi   % wSMI on windows
    end

    methods
        function [indicators, indicators_idx] = get_measures(obj,selected_measures)
            [ind, idx] = ismember(obj.all_names,selected_measures);
            all_indicators = {obj.cc,obj.corr,obj.coh,obj.im,obj.plv,obj.pli,obj.wpli,obj.mi,obj.wsmi};
            indicators = all_indicators(ind); indicators_idx = idx(ind);
        end

        % Associations computation
        function conn_mat = assoc_steps(obj,assoc_steps,indicators,x_bands,y_bands,x_hann,y_hann,fs,f_bands,select_info)
            % Common processing steps (association dependent)
            for step = assoc_steps
                switch step
                    case "cc";      obj = obj.cross_correlation(x_bands,y_bands,fs);
                    case "cpsd";    obj = obj.cross_spectrum(x_hann,y_hann,fs,f_bands);
                    case "cohy";    obj = obj.coherency;
                    case "info";    obj = obj.information(x_bands,y_bands,fs,select_info);
                end
            end

            % Specific processing steps
            conn_mat = zeros(length(indicators),length(f_bands));
            for ind = 1:length(indicators); conn_mat(ind,:) = indicators{ind}.function(obj); end
        end

        % PROCESSING STEPS BEFORE SPECIFIC ASSOCIATIONS COMPUTATION
        % Cross-correlation
        function obj = cross_correlation(obj,x,y,fs)
            n = size(x,2); maxlag_tp = floor(obj.cc.max_lag*fs);
            x_bar = mean(x,2); y_bar = mean(y,2); xy_std = std(x,0,2).*std(y,0,2);
            cross_correlogram_pos = zeros([maxlag_tp,size(x,1,3)]);
            cross_correlogram_neg = zeros([maxlag_tp,size(x,1,3)]);
            simple_corr = (1./(xy_std*n)).*sum((x-x_bar).*(y-y_bar),2);
            for lag = 1:maxlag_tp
                cross_correlogram_pos(lag,:,:) = (1./(xy_std*(n-lag))).*sum((x(:,1:end-lag,:)-x_bar).*(y(:,lag+1:end,:)-y_bar),2);
                cross_correlogram_neg(lag,:,:) = (1./(xy_std*(n-lag))).*sum((y(:,1:end-lag,:)-y_bar).*(x(:,lag+1:end,:)-x_bar),2);
            end
            obj.cross_corr = cat(1,flip(cross_correlogram_neg,1),permute(simple_corr,[2,1,3]),cross_correlogram_pos);
        end

        % Cross-spectrum
        function obj = cross_spectrum(obj,x,y,fs,f_bands)
            obj.f_bands = f_bands;
            [obj.cross_spec.sxy,obj.cross_spec.f] = cpsd(x,y,size(x,1),0,size(x,1),fs);
            [obj.cross_spec.sxx,~] = cpsd(x,x,size(x,1),0,size(x,1),fs);
            [obj.cross_spec.syy,~] = cpsd(y,y,size(x,1),0,size(x,1),fs);
        end

        % Coherency
        function obj = coherency(obj)
            coherency = obj.cross_spec.sxy./sqrt(obj.cross_spec.sxx.*obj.cross_spec.syy);
            obj.cohy = obj.mean_band(coherency);
        end

        % Mutual information and symbolic weighted version on windows
        function obj = information(obj,x,y,fs,select_info)
            miw = zeros([size(x,1) size(x,3)]); wsmiw = zeros(size(miw));
            x = (x - mean(x,2))./std(x,0,2); y = (y - mean(y,2))./std(y,0,2);
            for b = 1:size(x,1)
                tau_b = floor(obj.wsmi.tau(b)*fs);
                for w = 1:size(x,3)
                    xw = squeeze(x(b,:,w)); yw = squeeze(y(b,:,w));
                    if select_info(1); [~,miw(b,w)] = obj.mutual_information(xw,yw,obj.mi.bins,ones(obj.mi.bins)); end
                    if select_info(2)
                        x_hat = obj.symbolic_trnsfrm(xw,tau_b); y_hat = obj.symbolic_trnsfrm(yw,tau_b);
                        wsmiw(b,w) = obj.mutual_information(x_hat,y_hat,factorial(obj.wsmi.kernel),obj.wsmi.weights); 
                    end
                end
            end
            obj.wind_mi = miw; obj.wind_wsmi = wsmiw;
        end

        % Average of spectral content on narrow frequency bands of interest
        function bands = mean_band(obj,freq_spect)
            bands = zeros(size(obj.f_bands,1),size(freq_spect,2));
            for b = 1:size(obj.f_bands,1); bands(b,:) = mean(freq_spect(obj.cross_spec.f>=obj.f_bands(b,1) & obj.cross_spec.f<obj.f_bands(b,2),:),1); end
        end

        % Mutual information
        function [umi,bmi] = mutual_information(obj,x,y,bins,weights)
            edges = linspace(min([min(x),min(y)]), max([max(abs(x)),max(abs(y))+0.1]), bins+1);
            hxy = histcounts2(x,y,edges,edges); pxy = hxy/sum(hxy,"all");
            hx = histcounts(x,edges); px = hx/sum(hx);
            hy = histcounts(y,edges); py = hy/sum(hy);
            ex = 0; ey = 0; umi = 0;
            for i = 1:bins
                if px(i) > 0; ex = ex + px(i)*log10(1/px(i)); end
                if py(i) > 0; ey = ey + py(i)*log10(1/py(i)); end
                for j = 1:bins; if pxy(i,j) > 0; umi = umi + weights(i,j)*pxy(i,j)*log10(pxy(i,j)/(px(i)*py(j))); end; end
            end
            bmi = (2*umi) / (ex + ey);
        end

        % Symbolic transform
        function symbolic = symbolic_trnsfrm(obj,signal,tau)
            symbolic = zeros(size(signal,2)-tau*(obj.wsmi.kernel-1),1);
            for s = 1:length(symbolic)
                range = s:tau:s+(obj.wsmi.kernel-1)*tau;
                if      (signal(range(2))<signal(range(1))) && (signal(range(2))>signal(range(3))); symbolic(s) = 1;    % negative line
                elseif  (signal(range(3))<signal(range(1))) && (signal(range(3))>signal(range(2))); symbolic(s) = 2;    % U left
                elseif  (signal(range(1))<signal(range(2))) && (signal(range(1))>signal(range(3))); symbolic(s) = 3;    % bridge right
                elseif  (signal(range(1))<signal(range(3))) && (signal(range(1))>signal(range(2))); symbolic(s) = 4;    % U right
                elseif  (signal(range(2))<signal(range(3))) && (signal(range(2))>signal(range(1))); symbolic(s) = 5;    % positive line
                elseif  (signal(range(3))<signal(range(2))) && (signal(range(3))>signal(range(1))); symbolic(s) = 6;    % bridge left
                end
            end
        end

        % ASSOCIATION MEASURE FUNCTIONS
        % Cross correlation (absolute maximum of the cross-correlogram)
        function res = computeCC(obj); res = max(abs(squeeze(mean(obj.cross_corr,3))),[],1); end

        % Corrected cross-correlation
        function res = computeCORR(obj)
            maxlag_cc = (size(obj.cross_corr,1)-1)/2;
            corr_cc = zeros([maxlag_cc size(obj.cross_corr,2,3)]);
            for lag = 1:maxlag_cc; corr_cc(lag,:,:) = .5*(obj.cross_corr(lag+maxlag_cc+1,:,:)-obj.cross_corr(-lag+maxlag_cc+1,:,:)); end
            res = max(abs(squeeze(mean(corr_cc,3))),[],1);
        end

        % Coherence
        function res = computeCOH(obj); res = abs(mean(obj.cohy,2)); end

        % Imaginary coherence
        function res = computeIM(obj); res = abs(imag(mean(obj.cohy,2))); end

        % PLV
        function res = computePLV(obj)
            plv_freqs = abs(mean(exp(1i.*angle(obj.cross_spec.sxy)),2));
            res = obj.mean_band(plv_freqs);
        end

        % PLI
        function res = computePLI(obj)
            pli_freqs = abs(mean(sign(imag(obj.cross_spec.sxy)),2));
            res = obj.mean_band(pli_freqs);
        end

        % wPLI
        function res = computeWPLI(obj)
            wpli_freqs = abs(mean(imag(obj.cross_spec.sxy),2))./mean(abs(imag(obj.cross_spec.sxy)),2);
            res = obj.mean_band(wpli_freqs);
        end

        % Mutual information across windows
        function res = computeMI(obj); res = mean(obj.wind_mi,2); end

        % wSMI
        function res = computeWSMI(obj)
            bwsmi = (1/log10(factorial(obj.wsmi.kernel)))*obj.wind_wsmi;
            bwsmi(bwsmi<0) = 0; bwsmi(isnan(bwsmi)) = 0;
            res = mean(bwsmi,2);
        end
    end
end
