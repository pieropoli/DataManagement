clear
clc
close

% Code to parse data from obspyDMT to Matlab. It makes everything in a
% structure, one structure for each event with all the files.
% It performs rotation as 3c are read!
%
%
% Piero Poli. 1/6/16, MIT, Cambridge
% v 0.0.0

disp('add control for later phases!!!!')
disp('add save in station folder not event folder')
addpath ~/Seismic_Matlab_Functions/


%% parameter for the parsing

FB = [0.8 2]; % frequency band [hz]
newTau=0.01; % new frequqncy band for resampling the data [sec]
Win=200; % time saved after theoretical p wave arrival time [sec]
WinBeforeP=50; % time before P waaves used to save the data [sec]
type='vel';
preset=0;
inputDirectory= '/Volumes/MIT_01/DataPiero/Iquique2014/2014-01-04_2015-12-19'; % this is the directory in which you have the data
outDirectory='/Volumes/MIT_01/MAT/Iquique2014';
winSNR=20;
SNRlim=10;
distanceLim=[0,180]; % in deg, is the distance from eq to the station
Phase='P';
%% prepare catalog
mkdir(outDirectory)
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
for i1 =  1 : size(lon,2)
    
    %% get event information
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
    cnt = 1;
    
    
    if strcmp('vel',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH_VEL/*Z']);
    elseif strcmp('dis',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH/*Z']);
    elseif strcmp('raw',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{i1,1}) '/BH_RAW/*Z']);
    end
    datadir=dir(DIR);
   % datadir(1:2)=[];
    DIR = DIR(1:end-3);
    for i2 = 1  : size(datadir,1)
        
        station = char(datadir(i2, 1).name);
        z=char(station);z(end)='Z';
        e=char(station);e(end)='E';
        n=char(station);n(end)='N';
        if exist([char(DIR) '/' char(z)],'file') ==2 && exist([char(DIR) '/' char(e)],'file')  ==2  && exist([char(DIR) '/' char(n)],'file')  ==2
            
            
            sz=readsac([char(DIR) '/' char(z)]);
            se=readsac([char(DIR) '/' char(e)]);
            sn=readsac([char(DIR) '/' char(n)]);
            
            
            dateev=accuratetime;
            t0seisz=sz.jr*86400 + sz.hr*3600 + sz.mn*60 + sz.sec;
            t0seisn=sn.jr*86400 + sn.hr*3600 + sn.mn*60 + sn.sec;
            t0seise=se.jr*86400 + se.hr*3600 + se.mn*60 + se.sec;
            
            doyev = datevec2doy(dateev);
            t0ev=doyev*86400+ dateev(4)*3600 + dateev(5)*60 + dateev(6);
            dt0z= t0ev - t0seisz ;%;- sz.pointe; %% careful with this
            dt0e= t0ev - t0seise ;%- se.pointe; %% careful with this
            dt0n= t0ev - t0seisn ;%- sn.pointe; %% careful with this
            
            % get the distance from the earthqquake
            [sz.dist,sz.az]=distance(sz.elat,sz.elon,sz.slat,sz.slon);
            
            if sz.dist > min(distanceLim) && sz.dist <max(distanceLim) % distance control
                if 1/sz.tau*0.8*0.5 > max(FB) && 1/sn.tau*0.8*0.5 > max(FB) && 1/se.tau*0.8*0.5 > max(FB) % check the frequency range being valid
                    
                    
                    switch Phase
                        case 'P'
                            [tP,qP] = TravelTimeTaupPhasesDistance(sz.dist,'ttp+',sz.edep,'prem');
                            [tS,qS] = TravelTimeTaupPhasesDistance(sz.dist,'tts+',sz.edep,'prem');
                        case 'S'
                            [tP,qP] = TravelTimeTaupPhasesDistance(sz.dist,'tts+',sz.edep,'prem');
                    end
                    % get the right phase
                    if strcmp('P',Phase)==1;phn=1;else phn=find(strcmp(qP.ph,Phase)==1);end
                    if length(phn)>1;phn=phn(1);end
                    
                    %% check later phases
                    DT = tP-tP(phn);
                    F = find(DT>0&DT<Win);
                    laterphasename = qP.ph(F);
                    laterphasetime = tP(F) - tP(phn);
                    
                    %% isolate P wave window
                    t1z = floor((tP(phn)-WinBeforeP+dt0z)/sz.tau);
                    t1e = floor((tP(phn)-WinBeforeP+dt0e)/se.tau);
                    t1n = floor((tP(phn)-WinBeforeP+dt0n)/sn.tau);
                    
                    
                    t2 = t1z + floor((WinBeforeP + Win)/sz.tau);
                    
                    tnoise = t1z - floor(winSNR/sz.tau);
                    
                    if  t2<length(sz.trace) && t1z<length(sz.trace) && t1z>0 ...
                            && t1n<length(sn.trace) && t1n>0 ...
                            && t1e<length(se.trace) && t1e>0
                        
                        
                        if isempty(phn)==0
                            
                            tmpzbb = sz.trace(t1z:t2);
                            tmpnbb = sn.trace(t1n:t2);
                            tmpebb = se.trace(t1e:t2);
                            
                            % filtering
                            [a1,b1]=butter(2,FB*2*sz.tau);
                            sz.trace=filtfilt(a1,b1,sz.trace);
                            [a1,b1]=butter(2,FB*2*se.tau);
                            se.trace=filtfilt(a1,b1,se.trace);
                            [a1,b1]=butter(2,FB*2*sn.tau);
                            sn.trace=filtfilt(a1,b1,sn.trace);
                            
                            tmpz = sz.trace(t1z:t2);
                            tmpn = sn.trace(t1n:t2);
                            tmpe = se.trace(t1e:t2);
                            
                            sig=tmpz(floor(WinBeforeP/sz.tau):floor(WinBeforeP/sz.tau+winSNR/sz.tau));
                            noise = tmpz(1:floor(winSNR/sz.tau));
                            SNR = rms(abs(sig))/rms(abs(noise));                            %% HERE
                            if SNR>SNRlim
                                
                                %% interpolate tmp
                                t = 0:sz.tau:(length(tmpz)-1)*sz.tau;
                                t2 = 0:newTau:(length(tmpz)-1)*sz.tau;
                                outz = interp1(t,tmpz,t2);
                                outzbb = interp1(t,tmpzbb,t2);
                                t = 0:se.tau:(length(tmpe)-1)*se.tau;
                                t2 = 0:newTau:(length(tmpe)-1)*se.tau;
                                oute = interp1(t,tmpe,t2);
                                outebb = interp1(t,tmpebb,t2);
                                t = 0:sn.tau:(length(tmpn)-1)*sn.tau;
                                t2 = 0:newTau:(length(tmpn)-1)*sn.tau;
                                outn = interp1(t,tmpn,t2);
                                outnbb = interp1(t,tmpnbb,t2);
                                %% adjust the length ...
                                if length(outz)>(WinBeforeP+Win)/newTau
                                    N=length(outz)-(WinBeforeP+Win)/newTau;
                                    outz(end-N+1:end)=[];
                                    outzbb(end-N+1:end)=[];
                                elseif length(outz)<(WinBeforeP+Win)/newTau
                                    N = (WinBeforeP+Win)/newTau - length(outz);
                                    outz(end+N)=0;
                                    outzbb(end+N)=0;
                                end
                                
                                %% ...
                                if length(oute)>(WinBeforeP+Win)/newTau  
                                    N=length(oute)-(WinBeforeP+Win)/newTau;
                                    oute(end-N+1:end)=[];
                                    outebb(end-N+1:end)=[];
                                elseif length(oute)<(WinBeforeP+Win)/newTau 
                                    N = (WinBeforeP+Win)/newTau - length(oute);
                                    oute(end+N)=0;
                                    outebb(end+N)=0;
                                end
                                
                                %% ...
                                if length(outn)>(WinBeforeP+Win)/newTau 
                                    N=length(outn)-(WinBeforeP+Win)/newTau;
                                    outn(end-N+1:end)=[];
                                    outnbb(end-N+1:end)=[];
                                elseif length(outn)<(WinBeforeP+Win)/newTau 
                                    N = (WinBeforeP+Win)/newTau - length(outn);
                                    outn(end+N)=0;
                                    outnbb(end+N)=0;
                                end
                                
                                %% rotation...
                                [~,baz]=distance(sz.slat,sz.slon,sz.elat,sz.elon);
                                [R,T] = RotateComp(outn,oute,baz);
                                [Rbb,Tbb] = RotateComp(outnbb,outebb,baz);
                                %% save...
                                data.datazbb(cnt,:)=outzbb;
                                data.datatbb(cnt,:)=Tbb;
                                data.datarbb(cnt,:)=Rbb;
                                data.dataz(cnt,:) = outz;
                                data.datat(cnt,:) = T;
                                data.datar(cnt,:) = R;
                                data.snr(cnt) = SNR;
                                data.slat(cnt)=sz.slat;
                                data.slon(cnt)=sz.slon;
                                data.staname{cnt}=sz.staname;
                                data.comp{cnt}='3c';
                                data.evdist(cnt)=sz.dist;
                                data.evazimuth(cnt)=sz.az;
                                data.takeoff(cnt)=qP.takeoff(phn);
                                data.inciangle(cnt)=qP.inciangle(phn);
                                data.rayp(cnt)=qP.p(phn);
                                data.phase{cnt}=char(qP.ph(phn));
                                data.laterphases{cnt}=laterphasename;
                                data.laterphasetime{cnt}=laterphasetime;
                                data.arrivaltime(cnt)=tP(phn);
                                data.arrivaltimeS(cnt)=tS(1);
                                clear outz oute outn SNR s t t2 noise sig T R
                                cnt =cnt+1;
                                
                            end
                            
                        end
                    end
                    
                end % check the frequency range being valid
                
            end % distance control
            
            
        end
        
    end
    
    if isfield(data,'dataz')==1
    eval(['save ' char(outDirectory) '/data_' char(nameEvents{i1,1}) '.mat data'])
    disp(['save ' char(outDirectory) '/data_' char(nameEvents{i1,1}) '.mat data'])
    else
    disp('no data!')
    
    end
    clear data
    
    
end
