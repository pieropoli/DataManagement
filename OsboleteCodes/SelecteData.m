clear
clc

% This program read the output data from ParseData.m and create a single
% dataset for further analysis


DIR='/Users/pieropoli/MatData/GlobalM7+100km-20052015/data_P';
d=dir(DIR);

outdir = '/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/AlignedData/GlobalM7+100km-20052015';
mkdir(outdir)

for idir =  3 : length(d)
%% PARAMETERS

dd = dir([char(DIR) '/' char(d(idir, 1).name)]);
count=1;
for i1 = 3 : length(dd)
    
   load([char(DIR) '/' char(d(idir, 1).name) '/' char(dd(i1, 1).name)]) 
   data=s;
   Fp=1;

   if data.Prayinfo.takeoff(Fp)~=0 && data.dist<10 || data.dist>30

   Data2(count,1:length(data.trace))=data.trace;
   az(count)=data.az;
   takeoff(count)=data.Prayinfo.takeoff(Fp);
   dist(count)=data.dist;
   lo(count)=data.slon;
   la(count)=data.slat;
   snr(count)=s.SNR;
   count=count+1;
   end
end


out.tau=data.tau;
out.slon=lo;
out.slat=la;
out.magnitude=s.magnitude;
out.edep=data.edep;
out.elon=data.elon;
out.elat=data.elat;
out.dist=dist;
out.takeoff=takeoff;
out.Data=Data2; 
out.az=az;
out.beforePhase=s.WinBeforeP;
out.eventTime=s.eventtime;
out.ID=s.eventID;
out.snr=snr;

eval(['save ' char(outdir) '/Data_' char(s.eventID) '.mat out'])
clear out Data2 az  takeoff lo la dist snr

end