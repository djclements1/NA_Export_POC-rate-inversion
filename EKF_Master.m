%% This script is designed to first read in all data, and intialize the 
% Extended Kalmann Filter and Smoother. 
% Functions of each step are then called, and results are returned and
% saved in a .mat file. 

%% Step 1: read in all data files

% Need to make sure file directorries are data/filename

load('../DATA/Daily_avg_core_v2.mat')

load('DATA/NA_resp.mat')        % Incubation respiration
res = RESP;                % Change the name of data, to prevent duplicates. 
load('DATA/RESPIRE.mat')        % Respire trap data
load('DATA/NPP_struct.mat')     % 14C NPP
load('DATA/POC_obs.mat')        % LVISF pump observations
load ('DATA/NPP13_struct')      % 13C NPP
load('DATA/zooplankton_regrid') % Daily Zooplankton observations. 

%% Combine Resp from incubations and RESPIRE traps

lg_RESP = nanmean(cat(3,abs(RESP.gridded1_1),res.lg_RQ1_1),3);
sm_RESP = (abs(res.sm_griddedRQ1_1));


%% Initialize model boundaries
z_bnds(:,1) = [0,20,50,75,125,175,330];
z_bnds(:,2) = [20,50,75,125,175,330,500];
z_sctr = nanmean(z_bnds,2);
zwdt = z_bnds(:,2) - z_bnds(:,1); % width of depth bins
dt = 1; % time in days
T = 1:dt:30;
z = [20 50 75 125 175 330 500];
z_wdth = [20 30 25 50 50 155 170];
st = DATA.start'; % 6th of May
nd = DATA.end';   % 27th of May
ct = 0;

%% Label all variables
vars = {'Cs','Cs','Cs','Cs','Cs','Cs','Cs',...
       'Cl','Cl','Cl','Cl','Cl','Cl','Cl',...
       'wg','wg','wg','wg','wg','wg','wg',...
       'WL','WL','WL','WL','WL','WL','WL',...
       'J0','J0','J0','J0','J0','J0','J0',...
       'JL','JL','JL','JL','JL','JL','JL',...
       'B0','B0','B0','B0','B0','B0','B0',...
       'BL','BL','BL','BL','BL','BL','BL',...
       'B2P','B2P','B2P','B2P','B2P','B2P','B2P',...
       'BM2','BM2','BM2','BM2','BM2','BM2','BM2'};
       %'Bg','Be'};
var_longname = {'conc. small POC','conc. small POC','conc. small POC','conc. small POC','conc. small POC','conc. small POC','conc. small POC',...
       'conc. large POC','conc. large POC','conc. large POC','conc. large POC','conc. large POC','conc. large POC','conc. large POC',...
       'Cs sinking speed','Cs sinking speed','Cs sinking speed','Cs sinking speed','Cs sinking speed','Cs sinking speed','Cs sinking speed',...
       'Cl sinking speed','Cl sinking speed','Cl sinking speed','Cl sinking speed','Cl sinking speed','Cl sinking speed','Cl sinking speed',...
       'Cs production','Cs production','Cs production','Cs production','Cs production','Cs production','Cs production',...
       'Cl production','Cl production','Cl production','Cl production','Cl production','Cl production','Cl production',...
       'Cs respiration rate cst.','Cs respiration rate cst.','Cs respiration rate cst.','Cs respiration rate cst.','Cs respiration rate cst.','Cs respiration rate cst.','Cs respiration rate cst.',...
       'Cl respiration rate cst.','Cl respiration rate cst.','Cl respiration rate cst.','Cl respiration rate cst.','Cl respiration rate cst.','Cl respiration rate cst.','Cl respiration rate cst.',...
       'Aggregation rate cst.','Aggregation rate cst.','Aggregation rate cst.','Aggregation rate cst.','Aggregation rate cst.','Aggregation rate cst.','Aggregation rate cst.',...
       'Disaggregation rate cst.','Disaggregation rate cst.','Disaggregation rate cst.','Disaggregation rate cst.','Disaggregation rate cst.','Disaggregation rate cst.','Disaggregation rate cst.'...
       'Grazing rate cst.','Egestion rate cst.'};
       
units = {'mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3',...
       'mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3','mmol m^-3',...
       'm day^-1','m day^-1','m day^-1','m day^-1','m day^-1','m day^-1','m day^-1',...
       'm day^-1','m day^-1','m day^-1','m day^-1','m day^-1','m day^-1','m day^-1',...
       'mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1',...
       'mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1','mmol m^-3 day^-1',...
       'day^-1','day^-1','day^-1','day^-1','day^-1','day^-1','day^-1',...
       'day^-1','day^-1','day^-1','day^-1','day^-1','day^-1','day^-1',...
       'day^-1','day^-1','day^-1','day^-1','day^-1','day^-1','day^-1',...
       'm^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1','m^3 mmol^-1 day^-1'...
       'm^3 mmol^-1 day^-1','day^-1'};
      
%% Fill in gaps on initiation day - this is a linear interpolation between points

% Define our split in size 
OBS_small = DATA.micro' ;
OBS_large = DATA.small' + DATA.large' +  DATA.macro';


tpp = find(isnan(OBS_small(st,:)));
if~isempty(tpp)
    ct = ct+1;
    tt(ct) = st;
end

for idd = 1:length(tpp)
    if tpp(idd)==7
        OBS_small(st,tpp(idd)) = (0+OBS_small(st,tpp(idd)-1))/2;
        OBS_large(st,tpp(idd)) = (0+OBS_large(st,tpp(idd)-1))/2;
    elseif tpp(idd) ==1
        OBS_small(st,tpp(idd)) = (OBS_small(st,tpp(idd)+1)+OBS_small(tt(ct-1),tpp(idd)))/2;
        OBS_large(st,tpp(idd)) = (OBS_large(st,tpp(idd)+1)+OBS_large(tt(ct-1),tpp(idd)))/2;
    else
        OBS_small(st,tpp) = (OBS_small(st,tpp(idd)+1)+OBS_small(st,tpp(idd)-1))/2;
        OBS_large(st,tpp) = (OBS_large(st,tpp(idd)+1)+OBS_large(st,tpp(idd)-1))/2;
    end
end

OBS_large(OBS_large==0) = nan;
OBS_small(OBS_small==0) = nan;

%% Initialize all state terms
% Small POC - with obs
 x_int = OBS_small(st,:)';               % x is the total initital concentration
 p_int = (OBS_small(st,:)' .*0.5).^2;    % p is the total initial variance;
 OBS = OBS_small;
 ERR = (OBS_small.*0.5).^2;

% Large POC - with obs 
tpp = nanstd(OBS_large,[],1);
for ind = 1:length(z)
 OBS(:,end+1) = OBS_large(:,ind);
 ERR(:,end+1) = (OBS_large(:,ind).*0.5).^2;
 x_int(end+1,:) = OBS_large(6,ind);    % x is the total initital concentration
 p_int(end+1,:) = (OBS_large(6,ind).*0.5).^2;    % p is the total initial variance;
end

 wg = 1*ones(length(z),1);        % Set  Constant sinking speed 
 for ind = 1:length(z)
 OBS(:,end+1) =nan;                % No observations made
 ERR(:,end+1) = nan;
 x_int(end+1,:) = wg(ind);         % x is the total initital concentration
 p_int(end+1,:) = (wg(ind)).^2;    % 100% relative error;
 end

 WL = 8*ones(length(z),1);         % Constant Sinking speed
 for ind = 1:length(z)
 OBS(:,end+1) = nan;               % No observations made
 ERR(:,end+1) = nan;
 x_int(end+1,:) = WL(ind);         % x is the initial uniform sinking speed
 p_int(end+1,:) = (WL(ind)).^2;    % 100% relative error. 
 end

 
 J0 = NPP13C.alpha(:,6).*NPP14C.gridded(:,6);  % Split sm and lg by NPP13
 J0(3) = round((J0(2)+0)/2,2);                                  
 J0(isnan(J0)) = 10^-4;
 for ind = 1:length(z)
 OBS(:,end+1) = NPP13C.alpha(ind,:).*NPP14C.gridded(ind,:);
 ERR(:,end+1) = (NPP13C.alpha(ind,:).*NPP14C.gridded(ind,:) * 0.5).^2;
 x_int(end+1,:) = J0(ind);    % x is the total initital concentration
 p_int(end+1,:) = (J0(ind).*0.5).^2;    % p is the total initial error;
 end 

 JL = (1-NPP13C.alpha(:,6)).* NPP14C.gridded(:,6); 
 JL(3) = round((JL(2)+0)/2,2);
 JL(isnan(JL)) = 10^-4;
 for ind = 1:length(z)
 OBS(:,end+1) = (1-NPP13C.alpha(ind,:)).*NPP14C.gridded(ind,:);
 ERR(:,end+1) = (((1-NPP13C.alpha(ind,:)).*NPP14C.gridded(ind,:)) *0.5).^2;
 x_int(end+1,:) = JL(ind);    % x is the total initital concentration
 p_int(end+1,:) = (JL(ind).*0.5).^2;    % p is the total initial error;
 end
 
 sm_RESP(sm_RESP==0)=nan;
 B0 = sm_RESP(:,6);
 B0(isnan(B0)) = 0.05;
 for ind = 1:length(z)
 OBS(:,end+1) = sm_RESP(ind,:);
 ERR(:,end+1) = 0.1.^2;           % Absolut error;
 x_int(end+1,:) = B0(ind);        % x is the total initital concentration
 p_int(end+1,:) = 0.1.^2;         % p is the total initial error;
 end

 lg_RESP(lg_RESP==0) = nan;
 BL = lg_RESP(:,6);
 BL(isnan(BL)) = 0.02;
 for ind = 1:length(z)
 OBS(:,end+1) = lg_RESP(ind,:);
 ERR(:,end+1) = 0.1.^2;           % Absolute Error;
 x_int(end+1,:) = BL(ind);        % x is the total initital concentration
 p_int(end+1,:) = 0.1.^2;         % p is the total initial error;
 end
 
 B2P = zeros(length(z),1); % Aggregation
 B2P(B2P==0) = 0.003;       % Prior From Amaral et al. 2024
 for ind = 1:length(z)
 OBS(:,end+1) = nan;
 ERR(:,end+1) = nan;
 x_int(end+1,:) = B2P(ind);    % x is the total initital concentration
 p_int(end+1,:) = (B2P(ind).*0.5).^2;    % p is the total initial error;
 end
 
 BM2 = zeros(length(z),1);           % Disaggregation
 BM2(BM2==0) = 0.13;                 % Prior from Murnane 1996;
 for ind = 1:length(z)
 OBS(:,end+1) = nan;
 ERR(:,end+1) = nan;
 x_int(end+1,:) = BM2(ind);    % x is the total initital concentration
 p_int(end+1,:) = (BM2(ind).*0.5).^2;    % p is the total initial error;
 end

 %% DVM processing
 %Use Zp data
 %%% Step 1: all are perfectly known
 Zp = (zoo_grid.DVM*0.4)/12; %  Convert from Dry mass to Carbon Mass
 for ind = 1:31
     temp = find(Zp(ind,:)>0);
 Zg(ind) = zoo_grid.depths(min(temp));
 end

%%% Step 2: Variable Beta's
 Bg = 0.04;                % Rohr et al. 2022; 
 OBS(:,end+1) = nan;
 ERR(:,end+1) = nan;
 x_int(end+1,:) = Bg;    % x is the total initital concentration
 p_int(end+1,:) = (Bg).^2;    % p is the total initial error;
 vars{end+1} = 'Bg';

 % Determine the depth of grazing vs egestion!
 H = ones(7,1);
 cutoff = find(z==Zg(6));
 H(cutoff:end) = 0;
 Z = z_bnds(cutoff-1,2);
  
 % Integrate the total grazing
 INT_graze = (Bg*H'.*abs(Zp(6,:)).*OBS_small(6,:)) * z_wdth';
 
 % For just initialization - calculate the integrated eggestion as 30%
 % grazing

 INT_egg = 0.3*INT_graze;

% Calculate the intial specific egestion rate. 
 ZZ =  z_bnds(end,2) - Z;
 tot_ZOO = abs(Zp(6,:)).*(1-H') *z_wdth';
 Be = INT_egg/tot_ZOO; 
 OBS(:,end+1) = nan;
 ERR(:,end+1) = nan;
 x_int(end+1,:) = Be;    % x is the total initital concentration
 p_int(end+1,:) = (Be).^2;    % p is the total initial error;
 vars{end+1} = 'Be';

%% Initialize model errors
 q(1:7,1) = ones(7,1)*(1).^2; % POCs
 q(8:14,1) = ones(7,1)*(1).^2; %POCl
 q(15:21,1) = ones(7,1)*(1).^2; % Sinking Speed
 q(22:28,1) = ones(7,1)*(2).^2; % Sinking Speed

 tmp = [1 0.5 0.01 0.001 0.001 0.001 0.001];%flip(logspace(-3,0,7));
 q(29:35,1) = ones(7,1).*(tmp').^2;  % Small prod
 q(36:42,1) = ones(7,1).*(tmp').^2;  % Large prod

 q(43:49,1) = ones(7,1).*(0.1).^2;  % Small resp
 q(50:56,1) = ones(7,1).*(0.1).^2; % Large resp

 q(57:63,1) = ones(7,1)*(0.01).^2; % Aggregation
 q(64:70,1) = ones(7,1)*(0.05).^2; % Disagg

 q(71,1) = (0.01).^2; % Grazing beta  
 q(72,1) = (0.01).^2; % Eggestion Beta

 Q = diag(q);

%% Define mixed layer and Euphotic depth
 zmld = DA.mld03;
 zmld(zmld<0) = 0;
 zeu = DA.z1pc; % 1 pc Euphotic depth


 %% Kallman Filter
% This should be a function !!
[XM XP PM PP A0 W_et_s W_et_l] = EKF_POC(vars,z_bnds,zmld,Zp,Q,OBS,ERR,x_int, p_int);


 %% Smoother

 % This should also be a function!
[XN PN C] = RTSsmoother_POC(XM,XP,PM,PP,A0);

 %% Save model output
 model.day_of_month = T;
 model.depth_bin = z_bnds;
 model.depth_ctr = z_sctr;
 model.variable = vars;
 model.variable_longname = var_longname;
 model.units = units;
 model.data_obs = OBS;
 model.data_error = ERR;
 model.forward_estimate = XM;
 model.filter_estimate = XP;
 model.smoother_estimate = XN;
 model.forward_error = PM;
 model.filter_error= PP;
 model.smoother_error = PN;
 model.Small_entrain = W_et_s;
 model.Large_entrain = W_et_l;
 save('DATA/EKF_results','model')
