filename = '20150905_Mito.xlsx';
type = 'Mito'
[T,data,data_std,time_pts] = MetaParameters(filename,type)
%Which Cell lines to consider
cell_lines = {'WM164','A375','WM88','SKMEL28','WM793','SKMEL5','WM2664','WM983B','WM115','A2058'};
MetaPlots(T,data,data_std,time_pts,cell_lines,type)
