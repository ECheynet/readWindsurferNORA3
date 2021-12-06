function [data] = get2DSS(targetLat,targetLon,targetYear,targetMonth,targetDay,targetHour)
%
% [data] =
%  get2DSS(targetLat,targetLon,targetYear,targetMonth,targetDay,targetHour)
% reads the 2D wave elevation spectra from the WINDSURFER/NORA3 wave field hindcast with a  time
% resolution of 1 h and a spatial resolution of ca. 3 km. 
% The data set covers the period from 1993 to 2020. The wave hindcast are
% obtained with WAM wave model version cycle 4.7.0.
% The Norwegian Meteorological Institute provided the data freely.
% These are available on their Thredds server [1].
% The documentation is contained in the NetCDF files that are read
% through the OPeNDAP  framework.
% [1] https://thredds.met.no/thredds
%
% Input:
%    * targetLat: [1x1] or [1xNlat] double: target latitude where the data
%         should be extracted
%    * targetLon: [1x1] or [1xNlon] double: target longitude where the data
%         should be extracted
%    * targetYear: [1x1] or [1xNyear] double: year(s) from which the data are
%         taken
%    * targetMonth: [1x1] or [1xNmonth] double: month(s) from which the data are
%         taken
%    * targetDay: [1x1] or [1xNday] double: day(s) from which the data are
%         taken
%    * targetHour: [1x1] or [1xNhour] double: hour(s) from which the data are
%         taken
%
% Outputs:
%   * data: structure with the following fields
%       - time: [1x1] datetime
%     - 'lon': [1 x 1] double: longitude (deg)
%     - 'lat': [1 x 1] double: latitude (deg)
%     - 'f': [Nf x 1]: frequency vector from 0.03 Hz to 0.54 Hz
%     - 'theta': [Nt x 1]: direction vector with a resolution of 15 deg.
%     - 'S': [Nt x Nf]: 2D directional spectrum of wave elevation in m^2/Hz
% 
% Author: E. Cheynet - UiB, Norway - last modified: 06-12-2021


%% Check month and day number + transform date into string
if targetMonth<10,    myMonth = ['0',num2str(targetMonth)];else    myMonth = num2str(targetMonth);end
if targetDay<10,    myDay = ['0',num2str(targetDay)];else    myDay = num2str(targetDay);end
if targetHour<10,    myHour = ['0',num2str(targetHour)];else    myHour = num2str(targetHour);end
myYear = num2str(targetYear);
%% Preallocation and initalisation
data = struct('time',[],'S',[],'theta',[],'f',[]);

urldat= ['https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/',...
    '/',myYear,'/',myMonth,'/SPC',myYear,myMonth,myDay,'00.nc'];

time0 = seconds(ncread(urldat,'time')) + datetime([myYear,myMonth,myDay],'InputFormat','yyyyMMdd');
[~,indHour]=min(abs(time0-datetime([myYear,myMonth,myDay,targetHour],'InputFormat','yyyyMMddHH')));

lon00 = double(ncread(urldat,'longitude'));
lat00 = double(ncread(urldat,'latitude'));
offset = 0.5;
ind = find(lat00>=targetLat-offset &...
    lat00 <= targetLat+offset &...
    lon00 >=targetLon-offset &...
    lon00 <= targetLon+offset);
fprintf(['Number of locations used for interpolation: ',num2str(numel(ind)),' \n']);

freq = double(ncread(urldat,'freq'));
theta = double(ncread(urldat,'direction'));

N = numel(theta);
Nf = numel(freq);

%% read the spectra for each selected tile
S = zeros(N,Nf,numel(ind));
for ii=1:numel(ind)
S(:,:,ii) = ncread(urldat,'SPEC',[1,1,ind(ii),1,indHour],[inf inf, 1,1,1]); % direction,freq,lat,1,time
end

lon00 = lon00(ind);
lat00 = lat00(ind);

%% % Interpolate the data on the target location
data.S = zeros(N,Nf);
for ii=1:numel(theta)
    for jj=1:numel(freq)
        dummyVar = squeeze(double(S(ii,jj,:)));
        F_Var = scatteredInterpolant(lat00,lon00,dummyVar,'linear','none');
        data.S(ii,jj) = F_Var(targetLat,targetLon);
    end
end

%% store the output into the structure data
data.f = freq;
data.theta =theta;
data.lat = targetLat;
data.lon = targetLon;
data.time = time0(indHour);

%% add a value at 360 deg

if max(data.theta)<360
    data.theta(end+1) = 360;
    data.S(end+1,:) = mean(data.S([1 end],:));
    
    data.theta = [0;data.theta];
    data.S = [data.S(end,:);data.S];
end

% 


end