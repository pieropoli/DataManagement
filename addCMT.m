clear
clc
addpath /Users/pieropoli/Seismic_Matlab_Functions
globalcmt_update

%% load the data

DIR='invertedSpec2';
d=dir(DIR);

for idir =  3:length(d)

   % load data
   load([char(DIR) '/' char(d(idir, 1).name)]) 

   time= data.eventtime;
   [cmt,dist]=findcmt('time',time,'depth',data.evdep,'location',[data.evlat data.evlon],'magnitude',data.magnitude);
   out.cmt.cmt=cmt;
   out.cmt.dist=dist; % distance in cmt search
   out.scalarmoment=(cmt.scalarmoment*10^cmt.exponent)/10000000;
   eval(['save ' char(DIR) '/' char(d(idir, 1).name) ' out']) 
   dist
end
