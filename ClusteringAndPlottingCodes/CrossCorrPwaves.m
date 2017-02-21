clear
clc

% This program read the output data from ParseData.m and dot he correlaiton
% and shift of threaces.

winbefore=3; % before p wave
clim=0.5;
winCorr=4;  % window in which the correlation is evaluated



phase='P';
DIR='/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/DataManagement/test';
d=dir(DIR);

outdir = ['/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/AlignedData/PeruBrasil2015_TA/' char(phase)];
mkdir(outdir)

for idir = 3 : length(d)
%% PARAMETERS

dd = dir([char(DIR) '/' char(d(idir, 1).name)]);

for i1 = 1 : length(dd)
    
   load([char(DIR) '/' char(d(idir, 1).name) '/' char(dd(i1, 1).name)]) 
   data=s;
   if i1==3
       Data=zeros(length(dd)-2,(winCorr)/data.tau+winbefore/data.tau); 
       Data2=zeros(size(s.trace)); 
       az=zeros(1,length(dd)-2);
       tkoff=az;
       dist=az;
       lo=az;
       la=az;
   end

   data.trace=data.trace;
   if strcmp('P',phase)==1;
   Fp=1;
   else
   Fp=find(strcmp((phase),(data.Prayinfo.ph))==1);
   end
   Fp=Fp(1);
   
   tmp=data.trace(round(data.WinBeforeP-winbefore)/data.tau:round(data.WinBeforeP+winCorr)/data.tau);

   Data(i1-2,1:length(tmp))=tmp-(tmp(1)); 
   clear tmp

   tmp=data.trace;
   Data2(i1-2,1:length(tmp))=tmp.*tukeywin(length(tmp),.1)'; 
   Data2(i1-2,:)=Data2(i1-2,:)-mean(Data2(i1-2,:));
   
   az(i1-2)=data.az;
   tkoff(i1-2)=data.Prayinfo.takeoff(Fp);
   dist(i1-2)=data.dist;
   lo(i1-2)=data.slon;
   la(i1-2)=data.slat;

end


%% 1st alignment

[DataAlign,corrCoeff,delay] = MccAlignment(Data);
pol=mean(corrCoeff,1);
corrCoeff=abs(corrCoeff);
dt=mean(delay,1);
cc=mean(corrCoeff,1);
amp=rms(DataAlign,2);

F = find(abs(cc)<clim);
DataAlign(F,:)=[];
cc(F)=[];
Data2(F,:)=[];
dt(F)=[];
pol(F)=[];
az(F)=[];
dist(F)=[];
lo(F)=[];
la(F)=[];
tkoff(F)=[];
amp(F)=[];
clear F


if length(cc)>5
meanDT(1)=mean(abs(dt));
DataA=zeros(size(Data2));
DataNoNorm=DataA;
%% align
for i1 = 1 : size(Data2,1)
     
   out=delayTrace(Data2(i1,:),dt(i1)); 
   DataA(i1,:)=out./max(abs((DataAlign(i1,:))));%.*sign(pol(i1)); 
   DataAlign(i1,:)=DataAlign(i1,:).*sign(pol(i1));  
   DataNoNorm(i1,:)=out;  
    
end

%% 2nd alignment

clear corrCoeff delay
[DataAlign,corrCoeff,delay] = MccAlignment(DataAlign);
pol=mean(corrCoeff,1);
corrCoeff=abs(corrCoeff);
dt=mean(delay,1);
cc=mean(corrCoeff,1);


F = find(abs(cc)<clim);
cc(F)=[];
DataAlign(F,:)=[];
DataA(F,:)=[];
dt(F)=[];
pol(F)=[];
az(F)=[];
dist(F)=[];
lo(F)=[];
la(F)=[];
tkoff(F)=[];
amp(F)=[];
DataNoNorm(F,:)=[];
clear F


if length(cc)>5
meanDT(2)=mean(abs(dt));


for i1 = 1 : size(DataA,1)
    
    
   out=delayTrace(DataA(i1,:),dt(i1)); 
   DataA(i1,:)=out; 
    
    
    
end


%% 3rd alignment

[DataAlign,corrCoeff,delay] = MccAlignment(DataAlign);
pol=mean(corrCoeff,1);
corrCoeff=abs(corrCoeff);
dt=mean(delay,1);
cc=mean(corrCoeff,1);


F = find(abs(cc)<clim);
cc(F)=[];
DataAlign(F,:)=[];
dt(F)=[];
pol(F)=[];
az(F)=[];
dist(F)=[];
lo(F)=[];
la(F)=[];
tkoff(F)=[];
amp(F)=[];
DataA(F,:)=[];
DataNoNorm(F,:)=[];
clear F


if length(cc)>5
meanDT(3)=mean(abs(dt));

    
for i1 = 1 : size(DataA,1)
    
    
   out=delayTrace(DataA(i1,:),dt(i1)); 
   DataA(i1,:)=out; 
    
    
    
end




%% 4th alignment

[DataAlign,corrCoeff,delay] = MccAlignment(DataAlign);
pol=mean(corrCoeff,1);
corrCoeff=abs(corrCoeff);
dt=mean(delay,1);
cc=mean(corrCoeff,1);


F = find(abs(cc)<clim);
cc(F)=[];
DataAlign(F,:)=[];
dt(F)=[];
pol(F)=[];
az(F)=[];
dist(F)=[];
lo(F)=[];
la(F)=[];
tkoff(F)=[];
amp(F)=[];
DataA(F,:)=[];
DataNoNorm(F,:)=[];
clear F


if length(cc)>5
meanDT(4)=mean(abs(dt));

    
for i1 = 1 : size(DataA,1)
    
    
   out=delayTrace(DataA(i1,:),dt(i1)); 
   DataA(i1,:)=out; 
    
    
    
end




%% 5th alignment

[~,corrCoeff,delay] = MccAlignment(DataAlign);
pol=mean(corrCoeff,1);
corrCoeff=abs(corrCoeff);
dt=mean(delay,1);
cc=mean(corrCoeff,1);


F = find(abs(cc)<clim);
cc(F)=[];
dt(F)=[];
pol(F)=[];
az(F)=[];
dist(F)=[];
lo(F)=[];
la(F)=[];
tkoff(F)=[];
amp(F)=[];
DataA(F,:)=[];
DataNoNorm(F,:)=[];
clear F


if length(cc)>5
meanDT(5)=mean(abs(dt));

    
for i1 = 1 : size(DataA,1)
    
    
   out=delayTrace(DataA(i1,:),dt(i1)); 
   DataA(i1,:)=out; 
    
    
    
end

%% polarity correction
ref = mean(DataA);
ref=ref(s.WinBeforeP/s.tau-winbefore/s.tau:s.WinBeforeP/s.tau+winCorr/s.tau);

for ic = 1 : size(DataA,1)
    
   c=corrcoef(ref,s.WinBeforeP/s.tau-winbefore/s.tau:s.WinBeforeP/s.tau+winCorr/s.tau);
   pol(ic)=sign(c(1,2)); 
   DataA(ic,:)=DataA(ic,:).*pol(ic); 
    
    
end

%% save
clear out

out.meanDT=meanDT;
out.alignedData=DataA;
out.tau=data.tau;
out.slon=lo;
out.slat=la;
out.magnitude=s.magnitude;
out.edep=data.edep;
out.elon=data.elon;
out.elat=data.elat;
out.coccoceff=cc;
out.dist=dist;
out.corrWindow=winCorr;
out.takoff=tkoff;
out.dataNonormalization=DataNoNorm; 
out.amplitude=amp;
out.az=az;
out.beforePhase=s.WinBeforeP;
out.eventTime=s.eventtime;
out.ID=s.eventID;

eval(['save ' char(outdir) '/alignedData_' char(s.eventID) '.mat out'])
clear out DataA amp DataAlign lo la dist meanDT Data pol tkoff az DataNoNorm cc


end
end
end
end
end
end