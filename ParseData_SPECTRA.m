clear
clc
% close all
% this program parse the data downloaded with Obspy (ver. 2016) and make
% spectral analysis of P or S waves
% Piero Poli, MIT
addpath ../src
addpath ../EnergyEstimation/src_energy/
addpath  /Users/pieropoli/Seismic_Matlab_Functions/
warning off
%% parameter for the parsing
% minimum length of the signal in sec to accept the loaded sac file
FB = [0.01 20]; % frequency band [hz]
newTau=0.01; % new frequqncy band for resampling the data [sec]
type='vel';
preset=300;
inputDirectory= '/Volumes/DataChile_2/Data/Illapel/2015-09-01_2016-01-31'; % this is the directory in which you have the data
SNRlim=2;
distanceLim=[0,7]; % in deg, is the distance from eq to the station
Phase='P';
dlog=0.001; % log space frequency df
dbinlog=0.05; % log space frequency bin
tbp = 3.5; % time bandwidth fft
kspec=7; % number of tapers fft
WinBeforeP=1;

snrband = [.05 0.1 0.5 2 5 10 20];
maxfreq=max(snrband); % maximum frequency analysed

%% parse the catalog to retrieve the events

outDirectory=('/Users/pieropoli/MatData/IllapelSpectra');
fid1 = fopen(char([char(inputDirectory) '/EVENTS-INFO/catalog_table.txt']),'r');
[attributes]=textscan(fid1,'%f %f %f %f %s %f %s %s %s','headerlines',5);
fclose(fid1);

nameEvents=attributes(8);
magnitude=attributes(6);
eventTime=attributes(5);
dep=attributes(4);
prem = prem;


for i1 =  1 :  size(attributes{1},1)
    
    % get duration and Fc for the event
    
    MoEstimate =10^((3/2)*magnitude{1,1}(i1)+16.1)*1e-7;
    ds1=.1*1e6; % stress drop Pa
    
    % depth of the event and Vs
    ddev = dep{1,1}(i1);
    hchk = find(abs(prem.depth-ddev) == min(abs(prem.depth-ddev)));
    vs=prem.vs(hchk(1))*1000;
    Fc1 = (16/7*ds1/MoEstimate)^(1/3)*0.9*vs; % lowest frequency considered
    Fc2 = maxfreq; % lowest frequency considered
    
    if magnitude{1,1}(i1) > 6;Fc1=0.01;end
    
    % window length
    Win =ceil(1/(Fc1))*4;
    if Win<5;Win=5;end
    
    if magnitude{1,1}(i1) > 6;Win=120;end
    
    
    % 3 Component readsac
    
    if strcmp('vel',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH_VEL']);
    elseif strcmp('dis',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH']);
    elseif strcmp('raw',type)==1
        DIR = ([char(inputDirectory) '/' char(nameEvents{1,1}{i1,1}) '/BH_RAW']);
    end
    
    
    system(['ls ' char(DIR) '/**HZ >tmplist'])
    fid1 = fopen('tmplist','r');
    [stations]=textscan(fid1,'%s');
    stations=stations{1};
    datadir=dir(DIR);
    fclose(fid1)
    
    
    for i2 = 1  : length(stations)
        
        
        % check if the three components exists
        z=stations(i2);

            
            % load files
            sz=readsac(char(z));
            %             se=readsac(char(e));
            %             sn=readsac(char(n));
            
            
            % synchro
            t0seis=sz.jr*86400 + sz.hr*3600 + sz.mn*60 + sz.sec;
            tev = char(eventTime{1,1}(i1));
            dateev=obspyTimeToMat(tev);
            doyev = datevec2doy(dateev);
            t0ev=doyev*86400+str2num(tev(12:13))*3600 + str2num(tev(15:16))*60 + str2num(tev(18:19));
            dt0z= t0ev - t0seis - preset;
            
            %% check that the highest frequency is not smaller than fc2
            
            if 1/(sz.tau)/2*0.8 >= Fc2
                
                % get the distance from the earthqquake
                [s.dist,s.az]=distance(sz.elat,sz.elon,sz.slat,sz.slon);
                % get the name with time information
                s.eventtime = (eventTime{1,1}{i1,1});
                s.edep=sz.edep;
                % get waves arrival times
                [tP,qP] = TravelTimeTaupPhasesDistance(s.dist,'ttp+',s.edep,'/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/EGF/src/chile');
                [tS,qS] = TravelTimeTaupPhasesDistance(s.dist,'tts+',s.edep,'/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/EGF/src/chile');
                
                % get the right phase
                if strcmp('P',Phase)==1;phn=1;tP=tP(phn);else phn=find(strcmp(qP.ph,Phase)==1);end
                
                if strcmp('S',Phase)==1;phn=1;tP=tS(phn);phn=1;end
                
                if length(phn)>1;phn=phn(1);end
                
                
                
                %% isolate P wave window
                t1 = floor(tP(1)/sz.tau+preset./sz.tau-WinBeforeP/sz.tau+dt0z/sz.tau);
                
                t2 = t1 + floor(Win/sz.tau);
                
                tnoise = t1 - floor(Win/sz.tau);
                
                
                
                if isempty(phn)==0
                    
                    % check distance and takeoff
                    if s.dist>min(distanceLim)&&s.dist<max(distanceLim)
                        
                        % get the magnitude
                        s.magnitude=magnitude{1,1}(i1);
                        s.eventID=char(nameEvents{1,1}{i1,1});
                        % filtering
                        if max(FB)>(1/sz.tau)/2
                            [a1,b1]=butter(2,[min(FB) (1/sz.tau/2)-1]*2*sz.tau);
                            sz.trace=filtfilt(a1,b1,sz.trace);
                            %                                se.trace=filtfilt(a1,b1,se.trace);
                            %                                sn.trace=filtfilt(a1,b1,sn.trace);
                        else
                            [a1,b1]=butter(4,FB*2*sz.tau);
                            sz.trace=filtfilt(a1,b1,sz.trace);
                            %                                se.trace=filtfilt(a1,b1,se.trace);
                            %                                sn.trace=filtfilt(a1,b1,sn.trace);
                        end
                        
                        
                        % isolate the phase analyzed and get the SNR
                        
                        
                        
                        
                        if tnoise>0 && t2<length(sz.trace) && t1<length(sz.trace) && t1>0 && Win < (tS(1) - tP(1))
                            
                            
                            tmpz = sz.trace(t1:t2);
                            noisez = sz.trace(tnoise:t1);
                            
                            %% do the spectra
                            [specsigz,frs,df] = mtspec(tmpz,sz.tau,tbp,kspec,1);
                            [specnoisez,frn] = mtspec(noisez,sz.tau,tbp,kspec,1);
                            
                            
                            % signal to noise ratio
                            f1=(1/Win); % lowest frequency of the spectrum that is of our interest
                            [~,uf] = min(abs(f1-snrband));
                            
                            for isnr = uf : length(snrband) -1
                                
                                fchk = (frs>=snrband(isnr)&frs<=snrband(isnr+1));
                                sigi(isnr) = mean(sqrt(specsigz(fchk)));
                                noisei(isnr) = mean(sqrt(specnoisez(fchk)));
                                SNRspec(isnr) = (sigi(isnr)./noisei(isnr));
                                if  sigi(isnr)/noisei(isnr) > SNRlim
                                    SNRchek(isnr) = 1;
                                else
                                    SNRchek(isnr) = 0;
                                end
                            end
                            fg = find(sigi==0);sigi(fg)=[];noisei(fg)=[];SNRchek(fg)=[];SNRspec(fg)=[];
                            
                            clear noisei sigi
                            % get the high signal noise part of the
                            % spectrum
                            
                            fchk = (frs>=snrband(uf)&frs<=snrband(end));
                            frs=frs(fchk);
                            specsigz=specsigz(fchk);
                            specnoisez=specnoisez(fchk);
                            
                            
                            if sum(SNRchek)==length(SNRchek)
                                
                                flog = log10(min(frs)):dlog:log10(max(frs));
                                logspecz = interp1(log10(frs),log10(specsigz),flog);
                                logspecnoise = interp1(log10(frs),log10(specnoisez),flog);
                                
                                
                                % interpolate tmp
                                t = 0:sz.tau:(length(tmpz)-1)*sz.tau;
                                t2 = 0:newTau:(length(tmpz)-1)*sz.tau;
                                
                                outz = interp1(t,tmpz,t2);
                                clear t1 t2
                                s.tau=newTau;
                                s.WinBeforeP=WinBeforeP;
                                s.WinAfterP=Win;
                                s.z=outz;
                                clear outz
                                s.slon=sz.slon;
                                s.slat=sz.slat;
                                s.edep=sz.edep;
                                s.elon=sz.elon;
                                s.elat=sz.elat;
                                s.phase=Phase;
                                s.freqBandFilter=FB;
                                s.datecreated=date;
                                s.Prayinfo=qP;
                                s.Srayinfo=qS;
                                s.speclogbinz=logspecz;
                                s.speclogfreq=flog;
                                s.specFreq=frs;clear frs
                                s.Fcs=[Fc1 Fc2];
                                s.SNRspec=SNRspec;clear SNRspec
                                s.win=Win;
                                s.df=df;
                                %%  save the data
                                mkdir([char(outDirectory) '/' char(nameEvents{1,1}{i1,1})])
                                eval(['save ' char(outDirectory) '/' char(nameEvents{1,1}{i1,1}) '/' char(sz.staname) '_' num2str(sz.an) '_' num2str(sz.jr) '_phase_' char(Phase) '_' char(type) ' s'])
                                disp(['save ' char(outDirectory) '/' char(nameEvents{1,1}{i1,1}) '/' char(sz.staname) '_' num2str(sz.an) '_' num2str(sz.jr) ' dist=' num2str(s.dist) ' takeoff=' num2str(qP.takeoff(phn))])
                                
                            else
                                disp(['SNR too low index=' num2str(i2)])
                                
                            end
                            
                            %end
                            
                        else
                            
                            disp('data too short')
                            
                        end
                        
                    else
                        disp(['no data distance range dist=' num2str(s.dist)])
                        
                    end
                    
                end
                
                
            else
                
                disp('beyond Fc2')
                
            end
            
            
        
        clear logspecnoise logspecz  logfreq flog logspecnoise logfreqnoise specnoise
    end
    
end
