clear
clc
close

% Code to parse data from obspyDMT to Matlab. It makes everything in a
% structure, one structure for each event
%
%
% Piero Poli. 1/6/16, MIT, Cambridge
% v 0.0.0


addpath ../src


%% parameter for the parsing

FB = [0.5 5]; % frequency band [hz]
newTau=0.01; % new frequency band for resampling the data [sec]
Win=100; % time saved after theoretical p wave arrival time [sec]
WinBeforeP=30; % time before P waaves used to save the data [sec]
type='vel';
preset=300;
inputDirectory= '/Users/pieropoli/Data/deepEGF/2015-01-01_2017-01-01/'; % this is the directory in which you have the data
outDirectory='~/MatData/DeepEGF';
winSNR=10;
SNRlim=2;
distanceLim=[0,180]; % in deg, is the distance from eq to the station
Phase='P';

%% parse the catalog to retrieve the events
mkdir(outDirectory);

%% prepare catalog
te = textread([char(inputDirectory) '/EVENTS-INFO/catalog.txt'],'%[^\n]');
cnt = 1;
for ind = 19 : 11 : size(te,1)
    
    tmp =te{ind+7};
    n = find(isspace(tmp)==1)+1;
    lat(cnt)= str2num(tmp(n:end));
    
    tmp =te{ind+8};
    n = find(isspace(tmp)==1)+1;
    lon(cnt)= str2num(tmp(n:end));
    
    tmp =te{ind+6};
    n = find(isspace(tmp)==1)+1;
    dep(cnt)= str2num(tmp(n:end));
    
    tmp =te{ind+4};
    n = find(isspace(tmp)==1)+1;
    magnitude(cnt)= str2num(tmp(n:end));
    
    tmp = te{ind+3};
    n = find(isspace(tmp)==1, 1, 'last' );
    tmp(1:n) = [];
    year =str2num(tmp(1:4));
    month  = str2num(tmp(6:7));
    day  = str2num(tmp(9:10));
    hour  = str2num(tmp(12:13));
    minut  = str2num(tmp(15:16));
    sec  = str2num(tmp(18:26));
    date(cnt,:) =  [year month day hour minut sec];
    
    
    tmp = te{ind+2};
    n = find(isspace(tmp)==1, 1, 'first' );
    nameEvents{cnt,:} = deblankl(tmp(n:end));
    
    tmp = te{ind+1};
    n = find(isspace(tmp)==1, 1, 'last' );
    catalogtype{cnt,:} = deblankl(tmp(n:end));
    
    
    cnt  =cnt+1;
end



%% loop over events
for i1 = 1 : size(lon,2)
    
    %% get event information
    data.evdep = dep(i1);
    data.evlat= lat(i1);
    data.evlon = lon(i1);
    data.magnitude = magnitude(i1);
    data.eventtime = date(i1,:);
    data.nameevent = nameEvents{i1,1};
    data.tau=newTau;
    data.freqBand=FB;
    data.winbeforephase=WinBeforeP;
    data.winafterphase=Win;
    data.type=type;
    data.catalog=catalogtype{i1,1};
    accuratetime = date(i1,:);
    
    cnt = 1; % counter
    if strcmp('vel',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH_VEL']);
    elseif strcmp('dis',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH']);
    elseif strcmp('raw',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH_RAW']);
    end
    datadir=dir(DIR);
    datadir(1:2)=[];
    
    for i2 = 1  : size(datadir,1)
        
        s=readsac([char(DIR) '/' char(datadir(i2, 1).name)]);                
        t0seis=s.jr*86400 + s.hr*3600 + s.mn*60 + s.sec;
        dateev=accuratetime;
        doyev = datevec2doy(dateev);
        t0ev=doyev*86400+ dateev(4)*3600 + dateev(5)*60 + dateev(6);
        dt0z= t0ev - t0seis ;%- s.pointe; %% careful with this        
        % get the distance from the earthqquake
        [s.dist,s.az]=distance(s.elat,s.elon,s.slat,s.slon);
        
        if s.dist > min(distanceLim) && s.dist <max(distanceLim) % distance control
            if 1/s.tau*0.8*0.5 > max(FB) % check the frequency range being valid

                switch Phase
                    case 'P'
                        [tP,qP] = TravelTimeTaupPhasesDistance(s.dist,'ttp+',s.edep,'prem');                        
                    case 'S'
                        [tP,qP] = TravelTimeTaupPhasesDistance(s.dist,'tts+',s.edep,'prem');
                end
                % get the right phase
                if strcmp('P',Phase)==1;phn=1;elseif strcmp('S',Phase)==1;phn=1;else phn=find(strcmp(qP.ph,Phase)==1);end
                
                if length(phn)>1;phn=phn(1);end

                %% check later phases
                DT = tP-tP(phn);
                F = find(DT>0&DT<Win);
                laterphasename = qP.ph(F);
                laterphasetime = tP(F) - tP(phn); 

                %% isolate P wave window
                t1 = floor((tP(phn)-WinBeforeP+dt0z)/s.tau);
                
                t2 = t1 + floor((WinBeforeP + Win)/s.tau);
                                
                if  t2<length(s.trace) && t1<length(s.trace) && t1>0
                    
                    if isempty(phn)==0
                        
                        % filtering
                        [a1,b1]=butter(2,FB*2*s.tau);
                        sbb = s.trace;
                        s.trace=filtfilt(a1,b1,s.trace);
                        
                        tmp = s.trace(t1:t2);
                        tmpbb=sbb(t1:t2);
                        
                        sig=tmp(floor(WinBeforeP/s.tau):floor(WinBeforeP/s.tau+winSNR/s.tau));
                        noise = tmp(1:floor(winSNR/s.tau));
                        SNR = rms(abs(sig))/rms(abs(noise));
                        
                        if SNR>SNRlim

                            % interpolate tmp
                            t = 0:s.tau:(length(tmp)-1)*s.tau;
                            t2 = 0:newTau:(length(tmp)-1)*s.tau;                       
                            out = interp1(t,tmp,t2);
                            outbb=interp1(t,tmpbb,t2);
                            % adjust the length
                            if length(out)>(WinBeforeP+Win)/newTau
                               N=length(out)-(WinBeforeP+Win)/newTau;
                               out(end-N+1:end)=[]; 
                               outbb(end-N+1:end)=[];
                            elseif length(out)<(WinBeforeP+Win)/newTau 
                               N = (WinBeforeP+Win)/newTau - length(out);
                               out(end+N)=0;  
                               outbb(end+N)=0;
                            end
                            data.data(cnt,:) = out;
                            databroadband.data(cnt,:) = outbb;
                            data.snr(cnt) = SNR;
                            data.slat(cnt)=s.slat;
                            data.slon(cnt)=s.slon;
                            data.staname{cnt}=s.staname;
                            data.comp{cnt}=s.kcomp;
                            data.evdist(cnt)=s.dist;
                            data.evazimuth(cnt)=s.az;
                            data.takeoff(cnt)=qP.takeoff(phn);
                            data.inciangle(cnt)=qP.inciangle(phn);
                            data.rayp(cnt)=qP.p(phn);
                            data.phase{cnt}=char(qP.ph(phn));
                            data.laterphases{cnt}=laterphasename;
                            data.laterphasetime{cnt}=laterphasetime;
                            data.arrivaltime(cnt)=tP(phn);
                            clear out SNR s t t2 noise sig 
                            cnt =cnt+1;
                        end
                    end
                end

            end % check the frequency range being valid
            
        end % distance control

    end
    
    % save event
    if isfield(data,'data')==1
    eval(['save ' char(outDirectory) '/data_' char(nameEvents{i1,1}) '.mat data'])
    disp(['save ' char(outDirectory) '/data_' char(nameEvents{i1,1}) '.mat data'])
    end
    clear data
    
end







