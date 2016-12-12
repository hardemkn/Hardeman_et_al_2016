function [] = MetaPlots(T,data,data_std,time_pts,cell_lines,type,varargin)
if isempty(varargin)
    fontsz = 18;
else
    fontsz = varargin(1);
end
groups = T.Properties.RowNames;
if ~isempty(cell_lines)
    cnt = 1;
    for i = 1:length(cell_lines)
        if sum(strcmpi(cell_lines{i},groups)) == 0
            str = sprintf('%s was not in the data');
            warning(str)
        else
            idx(cnt) = find(strcmpi(cell_lines{i},groups));
            cnt = cnt +1
        end
    end
    groups = {groups{idx}};
    data = data(idx,:);
    data_std = data_std(idx,:);
    T = T(idx,:);
end

switch type
    case 'Mito'
        %Plot the seahorse time course
        figure;
        hold on
        col = hsv(length(groups));
        for i = 1:length(groups)
            plot(time_pts,data(i,:),'color',col(i,:),'linewidth',5);
        end
        for i = 1:length(groups)
            shadedErrorBar(time_pts,data(i,:),data_std(i,:),{'o','color',col(i,:),'linewidth',4},.5);
        end
        legend(groups);
        title('Mitochondria Stress Test')
        xlabel('Time')
        ylabel('OCR')
        set(gca,'fontsize',fontsz)

        %Plot each cell line's 7 parameters as a bar plot
        % figure()
        % hold on;
        col = [241,86,102;66,191,235;209,119,180;93,108,168;170,212,110;66,188,125;0,0,0]./255;
        for i = 1:length(groups)
            figure()
            hold on
            h(1) = bar(1,T.non_mito_resp(i));
            h(2) = bar(2,T.basal_resp(i));
            h(3) = bar(3,T.atp_pro(i));
            h(4) = bar(4,T.proton_leak(i));
            h(5) = bar(5,T.max_resp(i));
            h(6) = bar(6,T.spare_cap(i));
            h(7) = bar(7,T.bio_health_idx(i));
            for j = 1:7
                set(h(j),'FaceColor',col(j,:));
            end
            if T.non_mito_resp(i) > 0
                h(1) = line([1,1],[T.non_mito_resp(i),T.non_mito_resp(i)+T.non_mito_resp_err(i)]);
            else
                h(1) = line([1,1],[T.non_mito_resp(i),T.non_mito_resp(i)-T.non_mito_resp_err(i)]);
            end
            if T.basal_resp(i) > 0
                h(2) = line([2,2],[T.basal_resp(i),T.basal_resp(i)+T.basal_resp_err(i)]);
            else
                h(2) = line([2,2],[T.basal_resp(i),T.basal_resp(i)-T.basal_resp_err(i)]);
            end
            if T.atp_pro(i) > 0
                h(3) = line([3,3],[T.atp_pro(i),T.atp_pro(i)+T.atp_pro_err(i)]);
            else
                h(3) = line([3,3],[T.atp_pro(i),T.atp_pro(i)-T.atp_pro_err(i)]);
            end        
            if T.proton_leak(i) > 0
                h(4) = line([4,4],[T.proton_leak(i),T.proton_leak(i)+T.proton_leak_err(i)]);
            else
                h(4) = line([4,4],[T.proton_leak(i),T.proton_leak(i)-T.proton_leak_err(i)]);
            end
            if T.max_resp(i) > 0
                h(5) = line([5,5],[T.max_resp(i),T.max_resp(i)+T.max_resp_err(i)]);
            else
                h(5) = line([5,5],[T.max_resp(i),T.max_resp(i)-T.max_resp_err(i)]);
            end        
            if T.spare_cap(i) > 0
                h(6) = line([6,6],[T.spare_cap(i),T.spare_cap(i)+T.spare_cap_err(i)]);
            else
                h(6) = line([6,6],[T.spare_cap(i),T.spare_cap(i)-T.spare_cap_err(i)]);
            end
            if T.bio_health_idx(i) > 0
                h(7) = line([7,7],[T.bio_health_idx(i),T.bio_health_idx(i)+T.bio_health_idx_err(i)]);
            else
                h(7) = line([7,7],[T.bio_health_idx(i),T.bio_health_idx(i)-T.bio_health_idx_err(i)]);
            end
            for j = 1:7
                set(h(j),'Color',col(j,:),'LineWidth',4);
            end
            title(groups{i})
            set(gca,'fontsize',fontsz)
            xticklabel_rotate([1:7],45,{'NMR','BR','AP','PL','MR','SC','BHI'},'fontsize',fontsz)
        end
        
        figure()
        plot(1,1,'.')
        abrv = {'NMR','BR','AP','PL','MR','SC','BHI'};
        str = sprintf('%s=Non-Mitochondria Resparation\n%s=Basal Resparation\n%s=ATP Production\n%s=Proton Leak\n%s=Max Resparation\n%s=Spare Capacity\n%s=Bioenergetic Health Index',abrv{1},abrv{2},abrv{3},abrv{4},abrv{5},abrv{6},abrv{7});
        h = text(0,1,str);
        h.FontSize = fontsz;
        set(gca,'visible','off')
        
        for j = 1:7
            figure('Position',[680 588 1000 800])
            hold on
            h = bar(table2array(T(:,j)));
            set(h,'FaceColor',col(j,:))
            for i = 1:length(groups)
                if table2array(T(i,j)) > 0
                    line([i,i],[table2array(T(i,j)),table2array(T(i,j))+table2array(T(i,j+7))]);
                else
                    line([i,i],[table2array(T(i,j)),table2array(T(i,j))-table2array(T(i,j+7))]);
                end
            end
            title(abrv{j})
            %set(gca,'fontsize',fontsz,'XTick',[1:length(groups)],'XTickLabel',groups)
            xticklabel_rotate([1:length(groups)],45,groups,'fontsize',fontsz)
        end
        
    case 'Glyco'   
        figure;
        hold on
        col = hsv(length(groups));
        for i = 1:length(groups)
            plot(time_pts,data(i,:),'color',col(i,:),'linewidth',5);
        end
        for i = 1:length(groups)
            if strcmp(groups{i},'Background');
                continue;
            end
            shadedErrorBar(time_pts,data(i,:),data_std(i,:),{'o','color',col(i,:),'linewidth',4},.5);
        end
        legend(groups);
        title('Glycolytic Function')
        xlabel('Time')
        ylabel('ECAR')
        set(gca,'fontsize',fontsz)
        
        col = [246,81,105;71,187,225;170,212,112;61,190,121]./250;
        for i = 1:length(groups)
            figure()
            hold on
            h(1) = bar(1,T.non_glyc_acid(i));
            h(2) = bar(2,T.glyc(i));
            h(3) = bar(3,T.glyc_cap(i));
            h(4) = bar(4,T.glyc_res(i));
            for j = 1:4
                set(h(j),'FaceColor',col(j,:));
            end
            if T.non_glyc_acid(i) > 0
                h(1) = line([1,1],[T.non_glyc_acid(i),T.non_glyc_acid(i)+T.non_glyc_acid_err(i)]);
            else
                h(1) = line([1,1],[T.non_glyc_acid(i),T.non_glyc_acid(i)-T.non_glyc_acid_err(i)]);
            end
            if T.glyc(i) > 0
                h(2) = line([2,2],[T.glyc(i),T.glyc(i)+T.glyc_err(i)]);
            else
                h(2) = line([2,2],[T.glyc(i),T.glyc(i)-T.glyc_err(i)]);
            end
            if T.glyc_cap(i) > 0
                h(3) = line([3,3],[T.glyc_cap(i),T.glyc_cap(i)+T.glyc_cap_err(i)]);
            else
                h(3) = line([3,3],[T.glyc_cap(i),T.glyc_cap(i)-T.glyc_cap_err(i)]);
            end        
            if T.glyc_res(i) > 0
                h(4) = line([4,4],[T.glyc_res(i),T.glyc_res(i)+T.glyc_res_err(i)]);
            else
                h(4) = line([4,4],[T.glyc_res(i),T.glyc_res(i)-T.glyc_res_err(i)]);
            end
            for j = 1:4
                set(h(j),'Color',col(j,:),'LineWidth',4);
            end
            title(groups{i})
            set(gca,'fontsize',fontsz)
            xticklabel_rotate([1:4],45,{'NGA','G','GC','GR'},'fontsize',fontsz)
        end
        
        figure()
        plot(1,1,'.')
        abrv = {'NGA','G','GC','GR'};
        str = sprintf('%s=Non-Glycolytic Acidification\n%s=Glycolysis\n%s=Glycolytic Capacity\n%s=Glycolytic Reserve',abrv{1},abrv{2},abrv{3},abrv{4});
        h = text(0,1,str);
        h.FontSize = fontsz;
        set(gca,'visible','off')
        
        for j = 1:4
            figure('Position',[680 588 1000 800])
            hold on
            h = bar(table2array(T(:,j)));
            set(h,'FaceColor',col(j,:))
            for i = 1:length(groups)
                if table2array(T(i,j)) > 0
                    line([i,i],[table2array(T(i,j)),table2array(T(i,j))+table2array(T(i,j+4))]);
                else
                    line([i,i],[table2array(T(i,j)),table2array(T(i,j))-table2array(T(i,j+4))]);
                end
            end
            title(abrv{j})
            set(gca,'fontsize',fontsz)
            xticklabel_rotate([1:length(groups)],45,groups,'fontsize',fontsz)
        end
        
end
