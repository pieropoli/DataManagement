clear
clc



winCorr=1; % size of the window for the correlation


load('/Users/pieropoli/Autmoatic_Parsing_Downloading_Events/DataManagement/test/data_20160824_1.mat')
Data = normalizeMy(data.data);
%Data=Data(:,2500:3500);

win = 5;
win = win/data.tau;

T = 2000 : win/2 : 10000;

%% FK

pos_0(1)=mean(data.slat);
pos_0(2)=mean(data.slon);


for i=1:length(data.slon)
    [dist(i),angN(i)]= distance(pos_0(1),pos_0(2),data.slat(i), data.slon(i));
    dist(i)=deg2km(dist(i));
    pos(i,1)=dist(i)*sind(angN(i));
    pos(i,2)=dist(i)*cosd(angN(i));
end


param.psi=0;
psi_rad=pi*param.psi/180;

Npos(:,1) = pos(:,1)*cos(psi_rad) + pos(:,2)*sin(psi_rad) ;
Npos(:,2) = pos(:,2)*cos(psi_rad) - pos(:,1)*sin(psi_rad) ;

Nb =  size(Npos,1) ;

p=-.2:.005:.2;

ux=p;
uy=p;

FB = [.8 2];

for ic = 1 : length(T) -1
    tmp = Data(:,T(ic):T(ic+1));

    freq = [0:ceil(length(tmp)/2)-1]/length(tmp)/data.tau ;

    datafft = fft(tmp') ;
    datafft = datafft(1:ceil(length(tmp)/2),:);

    for iux = 1 : length(ux)

        for iuy = 1 : length(uy)
            k             = 2*pi * [ux(iux) ; uy(iuy)] * freq ;
            phase         = - Npos * k ;
            replica       = exp(1i*phase.') ./ norm(exp(1i*phase.')) ;
            fk(:,iux,iuy) = mean(datafft .* replica,2) ;
        end

    end

    clear F
    F=find(freq>FB(1)&freq<FB(2));
    out=fk(F,:,:);
    m = squeeze(mean(abs(out),1));
    [maax,in] = max(m(:));
    [ii1,ii2] = ind2sub(size(m),in);
    MAX(ic) = maax;
    ppx(ic) = p(ii1);
    ppy(ic) = p(ii2);
    imagesc(ux,uy,squeeze(mean(abs(out),1)))
    OUT(ic,:,:) = squeeze(mean(abs(out),1));clear out
    pause
    close
    
    
end