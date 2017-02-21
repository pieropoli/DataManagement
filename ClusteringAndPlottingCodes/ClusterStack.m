clear
clc
close


%% given the data from ParseDataMatrix.m make cross correlation in a given window, it also make the clustes

t0t = 5; % time before the p waves
t1t = 5; % time after the p waves
crosscorrcutoff=0.95;
distlim=3;


%% load data

load('data_20160531_1.mat')

%% integrate data
[a1,b1]=butter(2,data.freqBand*2*data.tau);
for id = 1 : size(data.data,1)

   data.data(id,:) = filtfilt(a1,b1,cumsum(data.data(id,:))); 

    
end

%% make timing


t0 = floor((data.winbeforephase -t0t )/data.tau);
t1 = floor((data.winbeforephase + t1t)/data.tau);
tP = floor(data.winbeforephase /data.tau);

%% cross correlation
datatmp  = data.data(:,t0:t1);

%datatmp=datatmp(1:200,:);

[coeff,delay] = CrossCorrDelay(datatmp);
delay = delay.*data.tau;

%% cluster analysis

K = 1.001 - coeff;			% create dissimilarity matrix
K = K - diag(diag(K));			% remove diagonal (required format)
Y = squareform(K);          % transform to "pdist" vector format
Z = linkage(Y,'complete');  % linkage

clust = cluster(Z,'cutoff',1-crosscorrcutoff,'criterion','distance');

uni_clust = find(unique(clust));




%% stack clusters
cnt = 1;
for i1 = 1 : length(uni_clust)
    
    C(i1,:) = rand(1,3);
    
    F  = find(uni_clust(i1)==clust);
    
 
    d = data.evdist(F)

    f2 = find(abs(median(d)-d)<distlim);
    if length(f2) > 5 
    s1(cnt,:) = mean(data.data(F(f2),:));
    az(cnt) = mean(data.evazimuth(F(f2)));
    cnt =cnt+1;
    
    
    hold all
    plot(data.slon(F(f2)),data.slat(F(f2)),'.')
    
    end
end



%% plot


for iz = 1 : length(az)
    
    
   hold all
   plot(az(iz)+normalizeMy(s1(iz,:)).*20)
    
    
    
end

