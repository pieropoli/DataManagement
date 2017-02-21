clear
clc
close

addpath /Users/pieropoli/m_map/

%% given the data from ParseDataMatrix.m make cross correlation in a given window, it also make the clustes

t0t = 5; % time before the p waves
t1t = 15; % time after the p waves
crosscorrcutoff=0.85;

%% load data

load('/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/DataManagement/chileSouthNetGLOBAL/data_20161225_1.mat')

%% make timing


t0 = floor((data.winbeforephase -t0t )/data.tau);
t1 = floor((data.winbeforephase + t1t)/data.tau);
tP = floor(data.winbeforephase /data.tau);

%% filter [GET D]
[a1,b1]=butter(2,1*[.01 .05]*data.tau);

for ifi = 1 : size(data.data,1)
    
    data.data(ifi,:)=filtfilt(a1,b1,(data.data(ifi,:)));

end


%% cross correlation
datatmp  = data.data(:,t0:t1);

[~,corrCoeff,delay] = MccAlignment(datatmp);
for id = 1 : size(corrCoeff,1)
    
   DT = delay(1,1) - mean(delay(id,:)) ;
    
    % align the data
    if DT > 0
    datatmp = zeros(size(data.data(id,:)));
    datatmp(DT+1:end)=data.data(id,1:end-DT);
    
    elseif DT < 0
    dt = abs(DT);    
    datatmp = zeros(size(data.data(id,:)));
    datatmp(1:end-dt)=data.data(id,dt+1:end);
    
    elseif DT==0 
    datatmp = data.data(id,:);
    end
    
   DataAlign(id,:) = datatmp; clear datatmp
    
    
end

%datatmp=datatmp(1:200,:);

[coeff,delay] = CrossCorrDelay(DataAlign);
% delay = delay.*data.tau;
% 
% 



%% cluster analysis

K = 1.001 - coeff;			% create dissimilarity matrix
K = K - diag(diag(K));			% remove diagonal (required format)
Y = squareform(K);          % transform to "pdist" vector format
Z = linkage(Y,'complete');  % linkage

clust = cluster(Z,'cutoff',1-crosscorrcutoff,'criterion','distance');

subplot(224)
H=dendrogram(Z,0,'colorthreshold',1-crosscorrcutoff);

uni_clust = find(unique(clust));




%% plot clusters
subplot(2,2,[1 3])

t = 0 : data.tau : (size(data.data,2)-1)*data.tau;
for i1 = 1 : length(uni_clust)
    
    hold on
    C(i1,:) = rand(1,3);
    
    F  = find(uni_clust(i1)==clust);

    for i2 = 1 : length(F)
        
    plot(t-data.winbeforephase,i1+(data.data(F(i2),:)./max(abs(data.data(F(i2),:)))),'-','Color',C(i1,:) );

    end    
end
hold on
plot([0 0],[0 i1+1],'-k')
plot([-t0t -t0t],[0 i1+1],'r--','LineWidth',2)
plot([t1t t1t],[0 i1+1],'r--','LineWidth',2)
box on


%% plot map

subplot(222)
m_proj('stereographic','lat',data.evlat,'long',data.evlon,'radius',90);
m_grid('linewi',2,'tickdir','out');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
hold on
m_plot(data.evlon,data.evlat,'ok','MarkerFaceColor','r','MarkerSize',15)



for i1 = 1 : length(uni_clust)
    

    
    F  = find(uni_clust(i1)==clust);
    
    m_plot(data.slon(F),data.slat(F),'.','Color',C(i1,:))
   
end


%% other maps

figure(2)

m_proj('miller','lat',82);
m_coast('color','k');
m_grid('linestyle','none','box','fancy','tickdir','out');
hold on
m_scatter(data.slon(1:size(delay,1)),data.slat(1:size(delay,1)),50,median(delay),'filled','MarkerEdgeColor','k')
colorbar
caxis([-max(max(abs(delay))) max(max(abs(delay)))])