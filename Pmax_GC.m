function [Pmax, Qmax, m, a] = Pmax_GC(con_names, con_freq, con_tran, tran_length)
%Performs the Garrett & Cummins (2005) maximum theoretical power calculation 
%for a transect using constituent data & bathymetry data. To account for 
%nonuniform flow, the transect is considered in small, discrete segments. 
%The constituent and bathymetry input data must be interpolated
%to the midpoints of these segments. 
%
%Requires: t_tide.m (Pawlowicz, 2002)
%________________________________________________________________________________
%INPUT
%  NOTE: Ns = number of transect segments
%        Nc = number of tidal constituents
%
%  con_names: Ncx2 array of 2-character constituent names
%            ex: ['Q1';'O1';'K1';'N2';'M2';'S2';'M4';'M6']  
%
%  con_freq: Ncx1 array of constituent frequencies [1/hr]
%            ex: [.0372; .0387; .0418; .0790; .0805; .0833; .1610; .2415]
%
%  con_tran: NsxNc matrix of constituent data at each transect segment midpoint
%     NOTE: Constituent data must be interpolated from a field of data points
%           onto small, discrete transect segments of equal length that form the 
%           full transect when connected end-to-end. Data should be interpolated 
%           to the midpoints of these segments. 
%
%     Matrix Columns:
%      1 ---------------- latitude [deg] *
%      2 ---------------- longitude [deg] * 
%      3 ---------------- depth [m]
%      4:(Nc+3) --------- water level amplitude [m] **
%      (Nc+4):(2Nc+3) --- velocity ellipse major axis [m/s] **
%      (2Nc+4):(3Nc+3) -- velocity phase [deg] ** 
%      (3Nc+4):(4Nc+3) -- velocity ellipse inclination CCW from East [deg] **
%
%        * Points must be listed in the order in which they occur on the
%          transect. The direction in which the transect increases as points 
%          are added does not matter. 
%       ** Field is repeated for each constituent in the order of con_names
%
%  tran_length: length of full transect [m]
%     
%_______________________________________________________________________________
%OUTPUT
%  Pmax: maximum theoretical power across the transect [MW]
%  Qmax: maximum flow rate due to dominant constituent [m^3/s]
%  m: multiplying factor [-]
%  a: dominant constituent amplitude [m]
%
%_______________________________________________________________________________
%Author Info:
% Written by Alexandra C. Muscalus and Kevin A. Haas at the 
% Georgia Institute of Technology, May 2017
% Contact: amuscalus@gatech.edu, khaas@gatech.edu
%_______________________________________________________________________________
%t_tide.m Reference:
%     Downloaded at https://www.eoas.ubc.ca/~rich/
%     Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal 
%       "Harmonic Analysis Including Error Estimates in MATLAB 
%        using t_tide", Computers and Geosciences, 2002.

%--------------------------------------------------------------------------
%                          1. Set G&C constants 
%--------------------------------------------------------------------------

rho=1020; %density of water kg/m^3
gamma=0.22; %G&C coefficient
g=9.81; %gravitational acceleration m/s^2

%--------------------------------------------------------------------------
%                     2. Read in constituent data
%--------------------------------------------------------------------------

%Check that input is consistent for all constituents
num_cons=(size(con_tran, 2)-3)/4; %number of tidal constituents
if ~(floor(num_cons)==num_cons) || num_cons==0
    disp('Error: input not formatted propertly.')
    return;
end

%Find transect locations with NaN values
nan_list=[];
for i=1:size(con_tran,1) %loop through transect segments
   if sum(isnan(con_tran(i,:)))>=1
       nan_list=[nan_list; [i, con_tran(i,1), con_tran(i,2)]];      
   end
end
if size(nan_list,1)>=1
    nan_list_c=num2cell(nan_list);
    nan_list_c=[{'Index'}, {'Latitude'}, {'Longitude'}; nan_list_c];
    str1=['Warning: NaN values at ', num2str(length(nan_list(:,1))), ' point(s):'];
    disp(str1) 
    disp(nan_list_c)
    disp('Flow rates at these points will be zero.')
end

%Retrieve parameters
lat=con_tran(:,1);
lon=con_tran(:,2);
h=con_tran(:,3); %depth
amp=con_tran(:,4:(num_cons+3)); %water level amplitudes 
v_maj=con_tran(:,(num_cons+4):(2*num_cons+3)); %major axis velocities 
v_ph=con_tran(:, (2*num_cons+4):(3*num_cons+3)); %velocity phases
v_inc=con_tran(:, (3*num_cons+4):(4*num_cons+3)); %inclination angles
num_points=length(lat);

%Change points with NaNs such that flow will be zero & t_predic will run
if ~isempty(nan_list)
    h(nan_list(:,1))=0;
    amp(nan_list(:,1),:)=.0001; %small but not zero
    v_maj(nan_list(:,1),:)=0;
    v_ph(nan_list(:,1),:)=1;
    v_inc(nan_list(:,1),:)=1;
end

%--------------------------------------------------------------------------
%              3. Determine transect endpoints & segment lengths
%--------------------------------------------------------------------------

lat1=lat(1); %transect endpoints
lon1=lon(1);
lat2=lat(end);
lon2=lon(end);

if lat1>lat2 
    latswap=lat1;
    lat1=lat2;
    lat2=latswap;
    lonswap=lon1;
    lon1=lon2;
    lon2=lonswap;
end

spacing=floor(tran_length/num_points); %meters

%--------------------------------------------------------------------------
%        4. Calculate amplitude & multiplying factor for each segment
%--------------------------------------------------------------------------

%Find dominant constituent by largest amplitude (mean of transect points)
[~, dom_ind]=max(nanmean(amp));
dom_amp=amp(:,dom_ind); %dom amp at each point

%G&C Multiplying Factor calculation
r=amp./repmat(dom_amp,1,num_cons); %ratio of non-dominant amps to dom amp
r(:, dom_ind)=0; %remove ratio of dominant constituent (equal to 1)
r(isnan(r))=0; %remove nans 
mult_fac=1+(9/16)*sum((r.^2),2); %calculation from G&C 2005

%Set bad values to NaN
if ~isempty(nan_list)
    mult_fac(nan_list(:,1))=NaN;
    dom_amp(nan_list(:,1))=NaN;
end


%--------------------------------------------------------------------------
%         5. Calculate constituent flow rates  at each segment 
%--------------------------------------------------------------------------

%Get velocity component perpendicular to transect
tran_angle=atand((lat2-lat1)/(lon2-lon1)); %Angle of transect CCW from East
perp_angle=tran_angle-90; %Angle perpendicular to transect (flow angle)
diff_angle=abs(v_inc(:,:)-perp_angle);
perp_v_maj=(v_maj(:,:).*cosd(diff_angle)); %velocity in line with flow 

%Flow rate, Q, at each segment (i) for each constituent
Q_i=perp_v_maj.*repmat(h,[1,num_cons])*spacing; %Q for each constituent 

%--------------------------------------------------------------------------
%        6. Generate time series for Qdom_i at each segment
%--------------------------------------------------------------------------

%ten minute intervals for 3 days 
t_ten_3day=datenum(2009,1,1,1,0,0):datenum(0,0,0,0,10,0):datenum(2009,1,4,1,0,0);

%Initialize pointwise time series variables (each discrete point has ts)
ts_Qdom_i=zeros(length(t_ten_3day), num_points);

%Loop through transect points and get Qdom time series at each point
for x = 1:size(Q_i,1) %Loop through transect points
    ts_Qdom_i(:,x)=t_predic(t_ten_3day, con_names(dom_ind,:), con_freq(dom_ind),....
        [Q_i(x,dom_ind)', 0, v_ph(x,dom_ind)', 0], 'latitude', mean([lat1, lat2]));
end

%--------------------------------------------------------------------------
%             7. Get Qmax (single value) and Qmax_i (discrete) 
%--------------------------------------------------------------------------

%Get time series of Qall (year) and Qdom (3 days) summed across full transect
ts_Qdom=nansum(ts_Qdom_i, 2);

%Get Qmax and tmax (maximum dominant constituent flow & time it occurs)
[Qmax, tmax]=max(abs(ts_Qdom)); %get max Qdom value and time

%Get discrete Qmax (Qmax_i) at tmax
Qmax_i=ts_Qdom_i(tmax,:); %get pointwise Qdom at tmax

%Make primary flow positive for discrete data
Qmax_i=Qmax_i.*(nansum(Qmax_i)/abs(nansum(Qmax_i)));


%--------------------------------------------------------------------------
%              8. Calculate weighted a, weighted m, and Pmax
%--------------------------------------------------------------------------

%Weight mult_fac and amp by Qdom to get single values
weight_fac=abs(Qmax_i)./sum(abs(Qmax_i));
m=nansum(weight_fac.*mult_fac'); %weighted mult factor (m)
a=nansum(weight_fac.*dom_amp'); %weighted amplitude (a)

%Pmax (G&C 2005)
Pmax=abs(gamma*rho*g*sum(Qmax_i)*m*a/(10^6));
end




