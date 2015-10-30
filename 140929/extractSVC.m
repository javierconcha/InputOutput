function [column1,column2,column3,column4] = extractSVC(filepath)
fid = fopen(filepath);
s = textscan(fid, '%s', 'delimiter', '\n');
% search for "data=":
idx1 = find(strcmp('data= ',s{1}),1);
% and read from s{1}(idx1+1:idx2-1) the values using textscan again ...
data = s{1}(idx1+1:size(s{1},1)-1);
fclose(fid);

for index = 1:size(data,1)
    SVCdata(index,:) = str2num(cell2mat(data(index)));
end

column1 = SVCdata(:,1);
column2 = SVCdata(:,2);
column3 = SVCdata(:,3);
column4 = SVCdata(:,4);
