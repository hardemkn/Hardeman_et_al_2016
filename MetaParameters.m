function [T,data,data_std,time_pts] = MetaParameters(filename,type, varargin)
    if isempty(varargin)
        writefile = 0;
    end
    switch type
        case 'Mito' %For calculating metabolic parameters for the the mitochondria stress test
            orig_data = readtable(filename);
            %Separate out the group data
            groups = unique(orig_data.Group);
            %Separate out the time points
            time_pts = unique(orig_data.Time);
            %Create a matrix with all the mean OCR values for each cell type at each
            %time point.  Row = cell line, Column = time
            for i = 1:length(groups)
                for j = 1:length(time_pts)
                    idx_group = find(strcmp(orig_data.Group,groups{i}));
                    idx_tim= find(orig_data.Time == time_pts(j));
                    idx_samp = intersect(idx_tim,idx_group);
                    num_Reps(i,j) = length(idx_samp);
                    
                    Q1 =  quantile(orig_data.OCR(idx_samp),.25);
                    Q3 = quantile(orig_data.OCR(idx_samp),.75);
                    IQR = Q3-Q1;
                    include = find(orig_data.OCR(idx_samp)>Q1-1.5*IQR & orig_data.OCR(idx_samp)<Q3+1.5*IQR);
                    
                    %include = find(abs(orig_data.OCR(idx_samp)-mean(orig_data.OCR(idx_samp)))>2*std(orig_data.OCR(idx_samp))==0);
                    data(i,j) = mean(orig_data.OCR(idx_samp(include)));
                    data_std(i,j) = std(orig_data.OCR(idx_samp(include)));
%                     if length(include) ~= length(idx_samp)
%                         sprintf('group = %i,time = %i',i,j)
%                     end
                end
            end
            %Calculate the 7 parameters of mitochondria function
            non_mito_resp = mean(data(:,10:12),2);
            basal_resp = mean(data(:,1:3),2)-non_mito_resp;
            atp_pro = basal_resp + non_mito_resp - mean(data(:,4:6),2);
            proton_leak = mean(data(:,4:6),2)-non_mito_resp;
            max_resp = mean(data(:,7:9),2)-non_mito_resp;
            spare_cap = max_resp - basal_resp;
            %See http://www.ncbi.nlm.nih.gov/pubmed/24895057 for equation
            a = 1;b=1;c=1;d=1;
            bio_health_idx = log10(abs(((spare_cap.^a).*(atp_pro.^b))./((proton_leak.^c).*(non_mito_resp.^d))));

            %Find how the error propogates in the mean
            %Here (columns = 4 distinct time points measured, rows = cell type)
            for i = 1:4
                err(:,i) = sqrt(data_std(:,(i-1)*3+1).^2 + data_std(:,(i-1)*3+2).^2 + data_std(:,i*3).^2)/3;
            end
            %Calculate the error of each of the 7 parameter above
            non_mito_resp_err = err(:,4);
            basal_resp_err = sqrt(non_mito_resp_err.^2+err(:,1).^2);
            atp_pro_err = sqrt(basal_resp_err.^2+non_mito_resp_err.^2+err(:,2).^2);
            proton_leak_err = sqrt(err(:,2).^2+non_mito_resp_err.^2);
            max_resp_err = sqrt(err(:,3).^2+non_mito_resp_err.^2);
            spare_cap_err = sqrt(max_resp_err.^2+basal_resp_err.^2);
            bio_health_idx_err = ((proton_leak.^c).*(non_mito_resp.^d))./((spare_cap.^a).*(atp_pro.^b)).*sqrt(...
                ((a*spare_cap.^(a-1)).^2).*(spare_cap_err).^2+...
                ((b*atp_pro.^(b-1)).^2).*(atp_pro_err).^2+...
                ((-c*proton_leak.^(c-1)).^2).*(proton_leak_err).^2+...
                ((-d*non_mito_resp.^(d-1)).^2).*(non_mito_resp_err).^2);

            T = array2table([non_mito_resp,basal_resp,atp_pro,proton_leak,max_resp,spare_cap,bio_health_idx,...
                non_mito_resp_err,basal_resp_err,atp_pro_err,proton_leak_err,max_resp_err,spare_cap_err,bio_health_idx_err]);
            T.Properties.VariableNames = {'non_mito_resp','basal_resp','atp_pro','proton_leak','max_resp','spare_cap','bio_health_idx',...
                'non_mito_resp_err','basal_resp_err','atp_pro_err','proton_leak_err','max_resp_err','spare_cap_err','bio_health_idx_err'};
            T.Properties.RowNames = groups;
            if writefile == 1
                if isempty(varargin{2})
                    idx = strfind(filename,'.');
                    writetable(T,[char(filename(1:idx-1)) '.csv'],'WriteRowNames',1)
                else
                    writetable(T,varargin{2},'WriteRowNames',1)
                end
            end
        case 'Glyco'
            orig_data = readtable(filename);
            groups = unique(orig_data.Group);
            time_pts = unique(orig_data.Time);
            for i = 1:length(groups)
                for j = 1:length(time_pts)
                    idx_group = find(strcmp(orig_data.Group,groups{i}));
                    idx_tim= find(orig_data.Time == time_pts(j));
                    idx_samp = intersect(idx_tim,idx_group);
                    
                    Q1 =  quantile(orig_data.ECAR(idx_samp),.25);
                    Q3 = quantile(orig_data.ECAR(idx_samp),.75);
                    IQR = Q3-Q1;
                    include = find(orig_data.ECAR(idx_samp)>=Q1-1.5*IQR & orig_data.ECAR(idx_samp)<=Q3+1.5*IQR);
                    num_Reps(i,j) = length(include);
                    
                    %include = find(abs(orig_data.ECAR(idx_samp)-mean(orig_data.ECAR(idx_samp)))>2*std(orig_data.ECAR(idx_samp))==0);
                    data(i,j) = mean(orig_data.ECAR(idx_samp(include)));
                    data_std(i,j) = std(orig_data.ECAR(idx_samp(include)));
                     if length(include) ~= length(idx_samp)
                         sprintf('group = %i,time = %i',i,j)
                     end
                end
                %%%%%
                %%%%%
                %data(i,:) = data(i,:) - data(i,1);
            end
            non_glyc_acid = mean(data(:,[10:12]),2);
            glyc = mean(data(:,4:6),2)-non_glyc_acid;
            glyc_cap = mean(data(:,7:9),2)-non_glyc_acid;
            glyc_res = glyc_cap-glyc;

            for i = 1:3
                if i == 1
                    err(:,i) = sqrt(sum(data_std(:,[1:3,9:12]).^2,2))/6;
                else
                    err(:,i) = sqrt(data_std(:,(i-1)*3+1).^2 + data_std(:,(i-1)*3+2).^2 + data_std(:,i*3).^2)/3;
                end
            end

            non_glyc_acid_err = err(:,1);
            glyc_err = sqrt(non_glyc_acid_err.^2+err(:,2).^2);
            glyc_cap_err = sqrt(non_glyc_acid_err.^2+err(:,3).^2);
            glyc_res_err = sqrt(glyc_cap_err.^2+glyc_err.^2);

            T = array2table([non_glyc_acid,glyc,glyc_cap,glyc_res,...
                non_glyc_acid_err,glyc_err,glyc_cap_err,glyc_res_err]);
            T.Properties.VariableNames = {'non_glyc_acid','glyc','glyc_cap','glyc_res',...
                'non_glyc_acid_err','glyc_err','glyc_cap_err','glyc_res_err'};
            T.Properties.RowNames = groups;
            if writefile == 1
                if isempty(varargin{2})
                    idx = strfind(filename,'.');
                    writetable(T,[char(filename(1:idx-1)) '.csv'],'WriteRowNames',1)
                else
                    writetable(T,varargin{2},'WriteRowNames',1)
                end
            end
    end
end