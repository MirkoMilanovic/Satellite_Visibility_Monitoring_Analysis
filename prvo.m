%          Creates a stereographic plot of GPS satellite
%          orbits from an almanac. The plot is as seen
%          from the position (phi,lambda). All orbits
%          with elevation angle lower than a given
%          cut-off angle (mask) are omitted.
%          Almanac files most easily can be downloaded
%          from your own (high-end) receiver or
%          http://www.ngs.noaa.gov/CORS/Data.thml
%
%          An additional plot is created showing number
%          of visible satellites and when they are visible

% definisanje stila i velicine koriscenih tekstualnih zapisa
set(0,'DefaultTextFontName','Times');
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultTextFontSize',16);

% prvo je neophodno ucitavanje navigacionog fajla, mask ugla, 
% phi i lambda koordinata date lokacije
rinexe('prvo.15n','aau.nav');
almanac = 'aau.nav';
mask = 10;
phi = [45 0 0];
lambda = [21 0 0];

%citanje efemerida
fide = fopen(almanac,'r');
Eph = fread(fide,inf,'double');
m = length(Eph);
eph = reshape(Eph,21,m/21);

% transformacija date lokacije (phi,lambda,h) u (X,Y,Z)
Phi = dms2rad(phi(1),phi(2),phi(3));
Phi = Phi*180/pi;
Lambda = dms2rad(lambda(1),lambda(2),lambda(3));
Lambda = Lambda*180/pi;
[M(1,1),M(2,1),M(3,1)] = frgeod(6378137,298.257222101,Phi,Lambda,0);

% Racunanje (azimuta, elevacije) za svaki satelit na svakih 15 min.
% Koristimo po jednu efemeridu za svaki satelit PRN.
% Graficki prikaz koristimo samo za vizuelno posmatranje.
[prns, ind] = unique(eph(1,:));
az = ones(32,96)*inf;
el = ones(32,96)*inf;

for sat = [ind]
    start_time = eph(18,sat);
    j = 0;
    i = eph(1,sat);
    for time = start_time:900:start_time+86400
        S = satpos(time,eph(:,sat));
        j = j+1;
        [azimuth,elevation,distance] = topocent(M,S-M);
        az(i,j) = azimuth;
        el(i,j) = elevation;
    end
end

%%%%%%%%% stereografska projekcija  %%%%%%%%%%%%%
figure(1);
% komanda polarhg iscrtava koordinatne linije, 
% dodajemo krugove za 30 i 60 stepeni
polarhg([30 60])
XX = zeros(32,40)*inf; % kreiramo praznu matricu za podatke koje prikazujemo
YY = XX;

hold on
for k = 1:32
    if az(k,1) == 0, continue, end
    AZ = az(k,:);
    EL = el(k,:);
    % odstranjivanje podataka ispod ucitanog minimalnog ugla
    AZ(find(EL <= mask)) = nan; 
    EL(find(EL <= mask)) = nan;
    % konverzija u polarne koordinate
    xx = (90-EL).*cos(AZ*pi/180);
    yy = (90-EL).*sin(AZ*pi/180);
    XX(k,1:length(xx)) = xx;
    YY(k,1:length(yy)) = yy;
end % end-k
% kod prve koordinate se tekst pomera vertikalno (rastuci ka gore),
% kod druge koordinate se tekst pomera horizontalno (rastuci ka desno)
text(135,-95,{['Skyplot za poziciju (\phi, \lambda) = (' ...
    num2str(round(Phi)) '\circ, '  num2str(round(Lambda)) '\circ)']})
text(115,-45,{['Mask ugao ' num2str(mask) '\circ' ]}) %120
text(-120,-120,['Svi PRN osim  ' num2str(setdiff(1:32,prns)) ])
plot(XX',YY','linewidth',2)
hold off

print -depsc2 easy111

break

% priprema za prikaz vidljivosti  %%%%%%%%%%%%%%%%%

% biramo rezoluciju od 5 minuta,
% 24 sata puta 12 = 288 sto postaje opseg za j
satsum = zeros(1,288);
visible = zeros(2*(size(prns,2)+1),288);

for sat = [ind]
    Az = [];
    El = [];
    i = eph(1,sat);
    for j = 1:288
        time = 300*(j-1);
        S = satpos(time,eph(:,sat));
        [az,el,d] = topocent(M,S-M);
        if el > mask
            Az = [Az az];
            El = [El el];
            satsum(j) = satsum(j)+1;
            visible(2*i,j) = 1;
        end
    end
end

figure(2);
set(gca,'Fontsize',16);
area(satsum)
set(gca,'XTick',1:71:288)
set(gca,'XTickLabel',{'0','6','12','18','24'})
xlabel('GPS Vreme [h]')
ylabel('# vidljivost satelita')
title(['Mask ugao ' num2str(mask) '\circ'])

print -depsc2 easy112

figure(3);
set(gca,'Fontsize',16);
imagesc(flipud(visible)); colormap(gray)
set(gca,'XTick',1:71:288)
set(gca,'XTickLabel',{'0','6','12','18','24'})
set(gca,'YTick',-3:16:(2*(size(prns,2)+1)))
set(gca,'YTickLabel',{'2','8','16','24','32'});
xlabel('GPS Vreme [h]')
ylabel('PRN')
title('Pune linije ukazuju na vidljive satelite')
colormap summer

print -depsc2 easy113
%%%%%%%%%%%%%%%%%%%%% kraj %%%%%%%%%%%%%
