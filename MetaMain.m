fontsz = 18;
filename = '20150905_Mito.xlsx';
type = 'Mito'
[T,data,data_std,time_pts] = MetaParameters(filename,type)
%Which Cell lines to consider
cell_lines = {'WM164','A375','WM88','SKMEL28','WM793','SKMEL5','WM2664','WM983B','WM115','A2058'};
%MetaPlots(T,data,data_std,time_pts,cell_lines,type)

groups = T.Properties.RowNames;
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
T = T(idx,:);
parm = T(:,1:7);

% 
% cell_lines = {'WM164-G','A375-G','WM88-G','SKMEL28-G','WM793-G','SKMEL5-G','WM2664-G','WM983B-G','WM115-G','A2058-G'};
% filename = '20150729_MitoGlyco.xlsx';
% type = 'Glyco'
% [T,data,data_std,time_pts] = MetaParameters(filename,type)
% MetaPlots(T,data,data_std,time_pts,cell_lines,type)

cell_lines = {'WM164','A375','WM88','SKMEL28','WM793','SKMEL5','WM2664','WM983B','WM115','A2058'};

filename = 'KNH_20150314_Glyco.xlsx';
type = 'Glyco'
[T,data,data_std,time_pts] = MetaParameters(filename,type)
%MetaPlots(T,data,data_std,time_pts,cell_lines,type)

groups = T.Properties.RowNames;
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
T = T(idx,:);
parm = [parm,T(:,1:4)]; parm = table2array(parm);

temp = readtable('Melanoma_IC50.xlsx');
ic_fifty = table2array(temp(:,2));
groups = cell_lines;

%Make all parameters dimentionless ratios of maximum
for i = 1:size(parm,2)
    norm_parm(:,i) = parm(:,i)./max(parm(:,i));
end
norm_ic_fifty = ic_fifty./max(ic_fifty);

%Run PCA on all of the components
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm)
figure
scatter(SCORE(:,1),norm_ic_fifty,'linewidth',5)
for i = 1:size(parm,1)
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty,SCORE(:,1)]);
str = sprintf('Corr=%.3f',temp(1,2));
text(.6,.9,str,'fontsize',fontsz)
title('PCA with all metabolic components')
axis([-2,2,0,1.25]);
set(gca,'fontsize',fontsz)

%Run PCA on just the mitochondria parameters
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,1:7))
figure
scatter(SCORE(:,1),norm_ic_fifty,'linewidth',5)
for i = 1:size(parm,1)
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
title('PCA with mitochondria components')
axis([-2,2,0,1.25]);
set(gca,'fontsize',fontsz)

%Run PCA on just the glycolosis parameters
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,8:11))
figure
scatter(SCORE(:,1),norm_ic_fifty,'linewidth',5)
for i = 1:size(parm,1)
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
title('PCA with all glycolosis parameters')
axis([-1,1,0,1.25]);
set(gca,'fontsize',fontsz)


abrv = {'NMR','BR','AP','PL','MR','SC','BHI','NGA','G','GC','GR'};
figure()
for i = 1:size(parm,2)  
    subplot(3,4,i)
    scatter(norm_parm(:,i),norm_ic_fifty,'linewidth',4)
    temp = corr([norm_ic_fifty,norm_parm(:,i)])
    c(i) = temp(1,2);
    text(.4,.7,sprintf('Correlation:%.3f',c(i)),'fontsize',14);
    title(char(abrv{i}))
    if i == 1
        xlabel('PCA1')
        ylabel('IC50')
    end
    text(norm_parm(1,i),norm_ic_fifty(1),'WM164','fontsize',16);
    text(norm_parm(10,i),norm_ic_fifty(10),'A2058','fontsize',16);
    set(gca,'fontsize',fontsz)
end

%Ignoring the A2058
figure()
for i = 1:size(parm,2)  
    subplot(3,4,i)
    scatter(norm_parm(1:end-1,i),norm_ic_fifty(1:end-1),'linewidth',4)
    temp = corr([norm_parm(1:end-1,i),norm_ic_fifty(1:end-1)])
    c(i) = temp(1,2);
    text(.2,.3,sprintf('Correlation:%.3f',c(i)),'fontsize',14);
    title(char(abrv{i}))
    if i == 1
        xlabel('PCA1')
        ylabel('IC50')
    end
    text(norm_parm(1,i),norm_ic_fifty(1),'WM164','fontsize',16);
    text(norm_parm(9,i),norm_ic_fifty(9),'WM115','fontsize',16);
    set(gca,'fontsize',fontsz)
end

figure()
i = 9
scatter(norm_parm(1:end-1,i),norm_ic_fifty(1:end-1),'linewidth',4)
temp = corr([norm_parm(1:end-1,i),norm_ic_fifty(1:end-1)])
c(i) = temp(1,2);
text(.5,.3,sprintf('Correlation:%.3f',c(i)),'fontsize',14);
title(char(abrv{i}))
if i == 1
    xlabel('PCA1')
    ylabel('IC50')
end
text(norm_parm(1,i),norm_ic_fifty(1),'WM164','fontsize',16);
text(norm_parm(9,i),norm_ic_fifty(9),'WM115','fontsize',16);
set(gca,'fontsize',fontsz)





[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[8,9,10]))
figure
scatter(norm_ic_fifty(1:end-1),SCORE(1:end-1,1),'linewidth',4)
for i = 1:size(parm,1)
    text(norm_ic_fifty(i),SCORE(i,1),groups{i},'fontsize',16);
end

corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)])

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[7,9]))
figure
scatter(norm_ic_fifty(1:end-1),SCORE(1:end-1,1),'linewidth',4)
for i = 1:size(parm,1)
    text(norm_ic_fifty(i),SCORE(i,1),groups{i},'fontsize',16);
end

corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)])

%Go Fishing
max_incl = 0;
min_incl = 0;
max_exl = 0;
min_exl = 0;
for i = 2:11
    V = nchoosek(1:11,i);
    for j = 1:size(V,1)
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,V(j,:)));
        temp = corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)]);
        temp2 = corr([norm_ic_fifty,SCORE(:,1)]);
        if temp(1,2) > max_incl
            max_incl = temp(1,2)
            store_1 = [];
            store_1 = V(j,:);
        elseif temp(1,2) < min_incl
            min_incl = temp(1,2)
            store_2 = [];
            store_2 = V(j,:);
        end
        if temp2(1,2) > max_exl
            max_exl = temp2(1,2)
            store_3 = [];
            store_3 = V(j,:);
        elseif temp2(1,2) < min_exl
            min_exl = temp2(1,2)
            store_4 = [];
            store_4 = V(j,:);
        end
    end
end

%Run PCA on the maximally correlated
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[store_1]))
figure
scatter(SCORE(1:end-1,1),norm_ic_fifty(1:end-1),'linewidth',5)
for i = 1:size(parm,1)-1
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)]);
str = sprintf('Corr=%.3f',temp(1,2));
text(.1,.25,str,'fontsize',fontsz)
title('PCA w/ G and GR excluding A2058')
set(gca,'fontsize',fontsz)

%Run PCA on the maximally correlated
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[store_1]))
figure
scatter(SCORE(:,1),norm_ic_fifty,'linewidth',5)
for i = 1:size(parm,1)
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty,SCORE]);
str = sprintf('Corr=%.3f',temp(1,2));
text(.1,.8,str,'fontsize',fontsz)
title('PCA w/ GR and GC including A2058')
set(gca,'fontsize',fontsz)




%Run PCA on the maximally correlated
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[store_3]))
figure
scatter(SCORE(:,1),norm_ic_fifty,'linewidth',5)
for i = 1:size(parm,1)
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty,SCORE(:,1)]);
str = sprintf('Corr=%.3f',temp(1,2));
text(.1,.5,str,'fontsize',fontsz)
title('PCA w/ PL, NGA, G, GC, GR')
set(gca,'fontsize',fontsz)


%Run PCA on the maximally correlated
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[10,11]))
figure
scatter(SCORE(1:end-1,1),norm_ic_fifty(1:end-1),'linewidth',5)
for i = 1:size(parm,1)-1
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)]);
str = sprintf('Corr=%.3f',temp(1,2));
text(-.3,.2,str,'fontsize',fontsz)
title('PCA w/ GR and GC exclude A2058')
set(gca,'fontsize',fontsz)

%Run PCA on the maximally correlated
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_parm(:,[9,10]))
figure
scatter(SCORE(1:end-1,1),norm_ic_fifty(1:end-1),'linewidth',5)
for i = 1:size(parm,1)-1
    text(SCORE(i,1),norm_ic_fifty(i),groups{i},'fontsize',16);
end
xlabel('PCA1')
ylabel('IC50')
temp = corr([norm_ic_fifty(1:end-1),SCORE(1:end-1,1)]);
str = sprintf('Corr=%.3f',temp(1,2));
text(.1,.25,str,'fontsize',fontsz)
title('PCA w/ G and GC exclude A2058')
set(gca,'fontsize',fontsz)



temp = []
for i = 2:11
    V = nchoosek(1:11,i);
    temp(i) = size(V,1);
end
