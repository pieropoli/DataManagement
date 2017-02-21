clear
clc
close all

addpath /Users/pieropoli/Seismic_Matlab_Functions

%% distance bin

bin = 3;
smooth = 1; % in sec
%% load data

load('/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/DataManagement/test/data_20160824_1.mat')

windowSize = round(smooth/data.tau);

b = (1/windowSize)*ones(1,windowSize);
a=1;

%% smooth
t = 0 : data.tau : (size(data.data,2)-1)*data.tau;
t = t -data.winbeforephase;

for is = 1 : length(data.slon)
ddata(is,:) = filter(b,a,data.data(is,:).^2);
norm = max(ddata(is,t>0&t<20));
ddata(is,:) =ddata(is,:) ./norm;

end

%% binning ...
D = 0 : bin : 90;
out = zeros(length(D),size(ddata,2));
N = zeros(size(D));

for i1 = 1 : length(D)
 
   F  = find(data.evdist>D(i1)&data.evdist<=(D(i1)+bin*2)); 
   
   if length(F)>1
   out(i1,:)  = sqrt(mean(ddata(F,:)));  
   norm = max(out(i1,t>0&t<20));    
   out(i1,:)=out(i1,:)./norm;
   N(i1)=length(F);
   %text(D(i1)+bin/2,-25,num2str(length(F)))
   end 
 
end


subplot(4,1,[1 2 3])
hold on
imagesc(D+bin/2,t,(out)'); 

for i1 = 1 : length(D)
    
    [tP,~] = TravelTimeTaupPhasesDistance(D(i1)+bin/2,'ttp+',data.evdep,'prem');
    [tS,~] = TravelTimeTaupPhasesDistance(D(i1)+bin/2,'tts+',data.evdep,'prem');
    
    plot(ones(length(tP)).*(D(i1)+bin/2),tP-tP(1),'ok','MarkerFaceColor','k')
    plot(ones(length(tS)).*(D(i1)+bin/2),tS-tP(1),'ok','MarkerFaceColor','g')
end


ylim([-30 max(t)-data.winbeforephase])
xlim([min(D) max(D)+bin])
caxis([0 .7])
colormap blue2red
ylabel('Time to P [s]')
title([char(data.eventtime) ' depth=' num2str(data.evdep) ' km M=' num2str(data.magnitude)])

subplot(4,1,4)
plot(D+bin/2,N);
xlim([min(D) max(D)+bin])
ylabel('# sta')
xlabel('Distance [deg]')