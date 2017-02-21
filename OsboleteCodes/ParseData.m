clear
clc
disp('modified the SNR! PLEASE CROOS CHECK')
% this program parse the data downloaded with Obspy (ver. may 2015)
% Piero Poli, MIT
addpath ../src
addpath /Users/pieropoli/Seismic_Matlab_Functions

warning off
%% parameter for the parsing
% minimum length of the signal in sec to accept the loaded sac file
FB = [0.01 1]; % frequency band [hz]
newTau=0.01; % new frequqncy band for resampling the data [sec]
Win=150; % time saved after theoretical p wave arrival time [sec]
WinBeforeP=20; % time before P waaves used to save the data [sec]
type='vel';
outtype='vel';
preset=300;
inputDirectory= '/Users/pieropoli/MatData/2016-04-16/2016-04-16_2016-04-17/'; % this is the directory in which you have the data
winSNR=10;
SNRlim=9;
distanceLim=[30,90]; % in deg, is the distance from eq to the station
Phase='P';

%% parse the catalog to retrieve the events

outDirectory=(['/Users/pieropoli/MatData/Ecuador_April2016/' char(Phase) '_' char(type)]);
fid1 = fopen(char([char(inputDirectory) '/EVENTS-INFO/catalog_table.txt']),'r');
[attributes]=textscan(fid1,'%f %f %f %f %s %f %s %s %s','headerlines',5);
fclose(fid1);

nameEvents=attributes(8);
magnitude=attributes(6);
eventTime=attributes(5);
dep=attributes(4);

for i1 = 1 : size(attributes{1},1)
    % readsac
    
    if strcmp('vel',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH_VEL']);
    elseif strcmp('dis',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH']);
    elseif strcmp('raw',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH_RAW']);
    end
    datadir=dir(DIR);
    
    mkdir([char(outDirectory) '/' char(nameEvents{1,1}{i1,1})])
    
    for i2 = 3  : size(datadir,1)
        
        
       s=readsac([char(DIR) '/' char(datadir(i2, 1).name)]);
        
        
       t0seis=s.jr*86400 + s.hr*3600 + s.mn*60 + s.sec;
       tev = char(eventTime{1,1}(i1));
       dateev=obspyTimeToMat(tev);
       doyev = datevec2doy(dateev);
       t0ev=doyev*86400+str2num(tev(12:13))*3600 + str2num(tev(15:16))*60 + str2num(tev(18:19));
       dt0z= t0ev - t0seis - preset;
        
        
        

        %% change unit
        if strcmp(outtype,'dis')==1 && strcmp(type,'vel')
            s.trace=cumsum(s.trace).*s.tau;
            disp('integration')
            
        end
        
        
        clear B F T
        % check the length of the trace
        
        
        % get the distance from the earthqquake
        [s.dist,s.az]=distance(s.elat,s.elon,s.slat,s.slon);
        % get the name with time information
        s.eventtime = (eventTime{1,1}{i1,1});
        
        % get waves arrival times
        [tP,qP] = TravelTimeTaupPhasesDistance(s.dist,'ttp+',s.edep,'prem');
        [tS,qS] = TravelTimeTaupPhasesDistance(s.dist,'tts+',s.edep,'prem');

        % get the right phase
        if strcmp('P',Phase)==1;phn=1;else phn=find(strcmp(qP.ph,Phase)==1);end
        if length(phn)>1;phn=phn(1);end

        %% isolate P wave window
        t1 = floor(tP(1)/s.tau+preset./s.tau-WinBeforeP/s.tau+dt0z/s.tau);
                                
        t2 = t1 + floor(Win/s.tau);
                                
        tnoise = t1 - floor(Win/s.tau);
        

        
        if tnoise>0 && t2<length(s.trace) && t1<length(s.trace) && t1>0
            
            if isempty(phn)==0
                
                % check distance and takeoff
                if s.dist>min(distanceLim)&&s.dist<max(distanceLim)
                    
                    % get the magnitude
                    s.magnitude=magnitude{1,1}(i1);
                    s.eventID=char(nameEvents{1,1}{i1,1});
                    % filtering
                    if max(FB)>(1/s.tau)/2
                        [a1,b1]=butter(3,[min(FB) (1/s.tau/2)-1]*2*s.tau);
                        s.trace=filtfilt(a1,b1,s.trace);
                        
                    else
                        [a1,b1]=butter(2,FB*2*s.tau);
                        s.trace=filtfilt(a1,b1,s.trace);
                    end

                    tmp = s.trace(t1:t2);
                    
                    
                    sig=max(abs(tmp(1:floor(WinBeforeP/s.tau+winSNR/s.tau))));
                    noise = s.trace(tnoise:t1);
                    SNR = sig/mean(abs(noise));

                    if SNR>SNRlim
                        
                        % interpolate tmp
                        t = 0:s.tau:(length(tmp)-1)*s.tau;
                        t2 = 0:newTau:(length(tmp)-1)*s.tau;
                        
                        out = interp1(t,tmp,t2);
                        
                        % update the data
                        s.tau=newTau;
                        s.WinBeforeP=WinBeforeP;
                        s.WinAfterP=Win;
                        s.SNR=SNR;
                        s.trace=out;
                        s.phase=Phase;
                        s.freqBandFilter=FB;
                        s.datecreated=date;
                        s.Prayinfo=qP;
                        s.Srayinfo=qS;
                        s.meanAmp=rms(s.trace);
                        
                        
                        %%  save the data
                        eval(['save ' char(outDirectory) '/' char(nameEvents{1,1}{i1,1}) '/' char(s.staname) '_' num2str(s.an) '_' num2str(s.jr) '_phase_' char(Phase) '_' char(outtype) ' s'])
                        disp(['save ' char(outDirectory) '/' char(nameEvents{1,1}{i1,1}) '/' char(s.staname) '_' num2str(s.an) '_' num2str(s.jr) ' dist=' num2str(s.dist) ' takeoff=' num2str(qP.takeoff(phn))])
                        
                    else
                        disp(['SNR too low index=' num2str(i2)])
                    end
                    
                    %end
                    
                else
                    
                    disp(['no data distance range dist=' num2str(s.dist)])
                    
                end
                
                
                
            end
        else
            disp('data too short')
        end
        
    end
    
    
end






