classdef network
    properties
        nb_features % number of features to extract
    end

    methods
        function obj = network(nb_features)
            obj.nb_features = nb_features;
        end

        % Features extraction from connectivity matrix
        function params = compute_features(obj,matrix,m)
            params = zeros([size(matrix,1,2) obj.nb_features]);
            for i = 1:size(matrix,1)
                for b = 1:size(matrix,2)
                    sym_mat = squeeze(matrix(i,b,:,:)) + squeeze(matrix(i,b,:,:))';
                    for ch = 1:size(matrix,3); sym_mat(ch,ch) = nan; end
                    if strcmp(m.indicators{i}.name,'wsmi'), sym_mat = weight_conversion(sym_mat,'normalize'); end
                    switch obj.nb_features
                        case 1; params(i,b,:) = mean(sym_mat,"all",'omitnan');
                        case 5; params(i,b,:) = obj.case5(sym_mat);
                        otherwise; disp('Invalid number of features');
                    end
                end
            end
        end

        % 5 network parameters: density, clustering coefficient, characteristic path length, assortativity, eigenvector centrality
        function params = case5(obj,sym_mat)
            params = zeros(1,obj.nb_features);
            params(1) = obj.density(sym_mat);
            sym_mat(isnan(sym_mat)) = 0; if all(sym_mat == 0); return; end
            params(2) = mean(clustering_coef_wu(sym_mat),'all');
            params(3) = charpath(distance_wei(1./sym_mat),0,0);
            params(4) = obj.assortativity(sym_mat);
            params(5) = mean(eigenvector_centrality_und(sym_mat),'all');
        end

        % Density
        function dens = density(obj,mat)
            tot_weight = sum(mat,"all",'omitnan');
            max_mat = max(mat,[],"all")*size(mat,1)*(size(mat,2)-1);
            dens = tot_weight./max_mat;
            if isnan(dens); dens = 0; end
        end

        % Assortativity
        function ass = assortativity(obj,mat)
            str = strengths_und(mat);
            [i,j] = find(triu(mat,1)>0); stri = str(i); strj = str(j); K = length(i);
            sum1 = 0; sum2 = 0; sum3 = 0;
            for l = 1:size(stri,2)
                sum1 = sum1 + mat(i(l),j(l))*stri(l)*strj(l);
                sum2 = sum2 + sqrt(mat(i(l),j(l)))*((stri(l)+strj(l))/2);
                sum3 = sum3 + mat(i(l),j(l))*((stri(l)^2+strj(l)^2)/2);
            end
            ass = ( K*sum1 - sum2^2 ) / ( K*sum3 - sum2^2 );
        end
    end
end
