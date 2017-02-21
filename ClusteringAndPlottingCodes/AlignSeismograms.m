clear
clc
close

load('/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/DataManagement/test/data_20160824_1.mat')


% select window
imagesc(normalizeMy(data.data));
[x,y] = ginput(2);
x = round(x);
tmp = data.data(:,min(x):max(x));

[DataAlign,corrCoeff,delay] = MccAlignment(tmp);