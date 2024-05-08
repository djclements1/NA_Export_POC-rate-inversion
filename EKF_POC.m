function  [XM XP PM PP A0 W_et_s W_et_l] = EKF_POC(vars,z_bnds,zmld,Zp,Q,OBS,ERR,x_int, p_int)

T = 1:30;
st = 6; 
nd = 27;

z_sctr = nanmean(z_bnds,2);
zwdt = z_bnds(:,2) - z_bnds(:,1); % width of depth bins
dt = 1; % time in days
T = 1:dt:30;
z = [20 50 75 125 175 330 500];
z_wdth = [20 30 25 50 50 155 170];

for ind = 2:length(zmld)
    dzmld(ind) = (zmld(ind) - zmld(ind-1));
end

H_et = dzmld;
H_et(H_et<=0) = 0;
H_et(H_et>0) = 1;

for t = 2:30
    for ii = 1:7
        H_tmp = z_bnds(ii,1) - zmld(t-1);
        if H_tmp>0
            H1(t,ii) = 1;
        else 
            H1(t,ii) = 0;
        end
        H_tmp = zmld(t) - z_bnds(ii,1); 
        if H_tmp>0
            H2(t,ii) = 1;
        else 
            H2(t,ii) = 0;
        end
        H_tmp = z_bnds(ii,2) - zmld(t-1);
        if H_tmp>0
            H3(t,ii) = 1;
        else 
            H3(t,ii) = 0;
        end

        H_tmp = zmld(t) - z_bnds(ii,2); 
        if H_tmp>0
            H4(t,ii) = 1;
        else 
            H4(t,ii) = 0;
        end
    end
end

H_top = H_et'.*H1.*H2;
H_bot = H_et'.*H3.*H4;
H_Layer = H_et'.*(H_top)+(H_bot);
H_Layer(isnan(H_Layer)) = 0;
W_et_l = nan(31,7);
W_et_s = nan(31,7);

for ind = 1:31
    temp = find(Zp(ind,:)>0);
    Zg(ind) = z_bnds(min(temp),2);
end

% Initialize output parameters
Z_hat = nan(max(T),length(x_int));
Z_obs = nan(max(T),length(x_int));
XM = nan(length(T),length(x_int)); 
PM = nan(length(x_int),length(x_int),length(T));
XP = nan(length(T),length(x_int)); 
PP = nan(length(x_int),length(x_int),length(T));


for t = st:nd
    % At time of initial observations use the initialization values.      
        if t ==st 
            XM(t,:) = x_int';
            PM(:,:,t) = diag(p_int);        
            XP(t,:) = x_int';
            PP(:,:,t) = diag(p_int);
        end
    
        if(t>st)   
            if ~isnan(mean(OBS(t,:),'omitnan'))  % This checks to make sure there is at least one observation at this time. 

                Z = squeeze(OBS(t,:,:))'; 
                H = ones(length(x_int),1);
                H = diag(H); 
                R_i = squeeze(ERR(t,:,:))';

                tpp = find(isnan(Z));
                Z(tpp) = [];
                H(tpp,:) =[]; 
                R_i(tpp) = [];
                R_i = diag(R_i);

                K = PM(:,:,t)*H'*(H*PM(:,:,t)*H' + R_i)^-1;  % Kalman gain
                XP(t,:) = XM(t,:)'+ (K*(Z-H*XM(t,:)'));      % Update X estimates
                PP(:,:,t) = PM(:,:,t) - K*H*PM(:,:,t);       % Update Covariance matrix
   
            else   % % Updated estimates are the priors if no data. 
                XP(t,:) = XM(t,:);
                PP(:,:,t) = PM(:,:,t);
            end
        end

        id_wg = find(~cellfun(@isempty,strfind((vars),'wg')));
        id_WL = find(~cellfun(@isempty,strfind((vars),'WL')));
        id_J0 = find(~cellfun(@isempty,strfind((vars),'J0')));
        id_JL = find(~cellfun(@isempty,strfind((vars),'JL')));
        id_B0 = find(~cellfun(@isempty,strfind((vars),'B0')));
        id_BL = find(~cellfun(@isempty,strfind((vars),'BL')));
        id_BM2 = find(~cellfun(@isempty,strfind((vars),'BM2')));
        id_B2P = find(~cellfun(@isempty,strfind((vars),'B2P')));
        id_Be = find(~cellfun(@isempty,strfind((vars),'Be')));
        id_Bg = find(~cellfun(@isempty,strfind((vars),'Bg')));

        H = ones(7,1);
        cutoff = find(z==Zg(t));
        H(cutoff:end) = 0;

    % This solves for the Jacobian of equations f
    for ii = 1: length(x_int) % These are the equations
        for jj = 1:length(vars) % these are the variables
            if strcmp(vars(ii),'Cs')
                if ii == jj
                    A0(ii,jj,t) = (1 - XP(t,id_wg(ii))/zwdt(ii)* dt - XP(t,id_B0(ii))*dt...
                        -(2*(XP(t,id_B2P(ii))*XP(t,ii)))) - H_bot(t,ii)/(2*zwdt(ii))*dzmld(t)...
                        - abs(Zp(t,ii)*XP(t,id_Bg)* H(ii)); 
                    if ii >1 && jj>1
                        A0(ii,jj-1,t) = XP(t,id_wg(ii-1))/zwdt(ii-1) * dt+ + H_top(t,ii)/(2*zwdt(ii))*dzmld(t);
                    end
                    if jj<6
                        A0(ii,jj+1,t) = - H_top(t,ii)/(2*zwdt(ii))*dzmld(t);
                        A0(ii,jj+2,t) =  H_bot(t,ii)/(2*zwdt(ii))*dzmld(t);
                    end
             
                elseif strcmp(vars(jj),'Cl') 
                    if ii==jj-7     
                        A0(ii,jj,t) = (XP(t,id_BM2(ii))*dt);           
                    end
                elseif strcmp(vars(jj),'wg')
                    if ii==jj-14 
                        A0(ii,jj,t) = -(XP(t,ii)/zwdt(ii));        
                        if ii >1 && jj>1
                            A0(ii,jj-1,t)  = XP(t,ii-1)/zwdt(ii-1) * dt;
                        end   
                    end

                elseif strcmp(vars(jj),'J0')
                    if ii==jj-28 
                        A0(ii,jj,t) = 1;
                    end
             
             
                elseif strcmp(vars(jj),'B0')
                    if ii==jj-42
                        A0(ii,jj,t) = -(XP(t,ii));
                    end

                elseif strcmp(vars(jj),'B2P')  
                    if ii==jj-56
                        A0(ii,jj,t) = -(XP(t,ii)).^2;
                    end               

             
                elseif strcmp(vars(jj),'BM2')
                    if ii==jj-63
                        A0(ii,jj,t) = (XP(t,ii+7)); %POCL
                    end
              
             
                elseif strcmp(vars(jj),'Bg')
                    if ismember(ii,1:7)
                        A0(ii,jj,t) = - XP(t,ii)*abs(Zp(t,ii)*(H(ii)));
                    end
                end
            end

            if strcmp(vars(ii),'Cl') % If its not 'J' - in the simple case it is C
                if ii == jj
                    A0(ii,jj,t) = (1 - XP(t,id_WL(ii-7))/zwdt(ii-7) * dt...
                        - XP(t,id_BL(ii-7))*dt-(XP(t,id_BM2(ii-7))*dt))...
                        - H_bot(t,ii-7)/(2*zwdt(ii-7))*dzmld(t);       
                
                    if ii >8 && jj>8
                        A0(ii,jj-1,t) = XP(t,id_WL(ii-7-1))/zwdt(ii-7-1) * dt...
                            + H_top(t,ii-7)/(2*zwdt(ii-7))*dzmld(t);
                    end
               
                    if jj<13
                        A0(ii,jj+1,t) = - H_top(t,ii-7)/(2*zwdt(ii-7))*dzmld(t);
                        A0(ii,jj+2,t) =  H_bot(t,ii-7)/(2*zwdt(ii-7))*dzmld(t);
                    end
             
                elseif strcmp(vars(jj),'Cs') 
                    if ii==jj+7
                        A0(ii,jj,t) = +2*(XP(t,id_B2P(ii-7))*dt)*XP(t,ii-7);        
                    end
                elseif strcmp(vars(jj),'WL')
                    if ii==jj-14 
                        A0(ii,jj,t) = -(XP(t,ii)/zwdt(ii-7));        
                        if ii >8 && jj>8
                            A0(ii,jj-1,t)  = XP(t,ii-1)/zwdt(ii-8) * dt;
                        end
                    end
             
                elseif strcmp(vars(jj),'JL')
                    if ii==jj-28
                        A0(ii,jj,t) = 1;
                    end
             
                elseif strcmp(vars(jj),'BL')
                    if ii==jj-42
                        A0(ii,jj,t) = -XP(t,ii);
                    end

                elseif strcmp(vars(jj),'B2P')
                    if ii==jj-49
                        A0(ii,jj,t) = (XP(t,ii-7)).^2;
                    end                   
             
                elseif strcmp(vars(jj),'BM2')  
                    if ii==jj-56
                        A0(ii,jj,t) = -(XP(t,ii)); 
                    end
              
                elseif strcmp(vars(jj),'Be')
                    if ismember(ii,8:14)
                        A0(ii,jj,t) = +abs(Zp(t,ii-7)*(1-H(ii-7))); 
                    end
                end
            end   
            
            if strcmp(vars(ii),'wg')
                if ii == jj
                    A0(ii,jj,t) =1;
                end
            end
          
            if strcmp(vars(ii),'WL')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end

            if strcmp(vars(ii),'J0')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end

            if strcmp(vars(ii),'JL')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end
          
            if strcmp(vars(ii),'B0')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end
         
            if strcmp(vars(ii),'BL')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end
          
            if strcmp(vars(ii),'BM2')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end          
            if strcmp(vars(ii),'B2P')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end  

            if strcmp(vars(ii),'Bg')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end     
          
            if strcmp(vars(ii),'Be')
                if ii == jj
                    A0(ii,jj,t) = 1;
                end
            end     
        end
    end
%%% Determine the model evolution equations - forward model

% Set a heavy-side function for grazing H==1 if above zg

for ii = 1:7
    if ii == 1 
            
        Etop = 1/(z_bnds(ii,2) - z_bnds(ii,1))*dzmld(t)*H_top(t,ii)...
            *((0-XP(t,ii+1))/2); 
        Ebot = 1/(z_bnds(ii,2) - z_bnds(ii,1))*dzmld(t)*H_bot(t,ii)...
            *((XP(t,ii+2)-XP(t,ii))/2);
        W_et_s(t,ii) = Etop+Ebot;

        f(ii,t) = XP(t,ii)*( 1-  XP(t,id_wg(ii))/zwdt(ii) * dt...
            - XP(t,id_B0(ii))*dt-((XP(t,id_B2P(ii))*dt)*XP(t,ii)))...
            + W_et_s(t,ii) + XP(t,id_J0(ii))*dt + XP(t,ii+7)*XP(t,id_BM2(ii))...
            - XP(t,id_Bg)* XP(t,ii) * abs(Zp(t,ii)) * H(ii);
        
    elseif ii>1

        Etop = 1/(z_bnds(ii,2) - z_bnds(ii,1))*dzmld(t)*H_top(t,ii)...
            *((XP(t,ii-1)-XP(t,ii+1))/2);
        Ebot = 1/(z_bnds(ii,2) - z_bnds(ii,1))*dzmld(t)*H_bot(t,ii)...
            *((XP(t,ii+2)-XP(t,ii))/2);
        W_et_s(t,ii) = Etop+Ebot;

        f(ii,t) =XP(t,ii)*( 1-  XP(t,id_wg(ii))/zwdt(ii) * dt...
            - XP(t,id_B0(ii))*dt-((XP(t,id_B2P(ii))*dt)*XP(t,ii)))...
            + XP(t,id_J0(ii))*dt + XP(t,ii+7)*XP(t,id_BM2(ii))...
            + W_et_s(t,ii) +  XP(t,id_wg(ii-1))/zwdt(ii-1) * dt * XP(t,ii-1)...
            - XP(t,id_Bg) * XP(t,ii) * abs(Zp(t,ii))* H(ii);
    end
end

for ii = 8:14
    if ii == 8   
        Etop = 1/(z_bnds(ii-7,2) - z_bnds(ii-7,1))*dzmld(t)*H_top(t,ii-7)...
            *((0-XP(t,ii+1))/2);
        Ebot = 1/(z_bnds(ii-7,2) - z_bnds(ii-7,1))*dzmld(t)*H_bot(t,ii-7)...
            *((XP(t,ii+2)-XP(t,ii))/2); 
        W_et_l(t,ii-7) = Etop+Ebot;

        f(ii,t) = XP(t,ii)*(1 - XP(t,id_WL(ii-7))/zwdt(ii-7)*dt ...
            - XP(t,id_BL(ii-7))*dt- XP(t,id_BM2(ii-7))*dt) ...
            + W_et_l(t,ii-7) + XP(t,id_B2P(ii-7))*dt*XP(t,ii-7).^2 + XP(t,id_JL(ii-7))...
            + XP(t,id_Be) * abs(Zp(t,ii-7))*(1-H(ii-7));
    elseif ii>8
        Etop = 1/(z_bnds(ii-7,2) - z_bnds(ii-7,1))*dzmld(t)*H_top(t,ii-7)...
            *((XP(t,ii-1)-XP(t,ii+1))/2);
        Ebot = 1/(z_bnds(ii-7,2) - z_bnds(ii-7,1))*dzmld(t)*H_bot(t,ii-7)...
            *((XP(t,ii+2)-XP(t,ii))/2);
        W_et_l(t,ii-7) = Etop+Ebot;

            
        f(ii,t) = XP(t,ii)* (1 - XP(t,id_WL(ii-7))/zwdt(ii-7)*dt ...
            - XP(t,id_BL(ii-7))*dt- XP(t,id_BM2(ii-7))) ...
            + XP(t,id_B2P(ii-7))*dt*XP(t,ii-7).^2 + XP(t,id_JL(ii-7))...
            + W_et_l(t,ii-7) + XP(t,id_WL(ii-7-1))/zwdt(ii-7-1)*dt * XP(t,ii-1)...
            + XP(t,id_Be) * Zp(t,ii-7)*(1-H(ii-7));
    end
end
    
for ii = 15:72
    f(ii,t) =XP(t,ii);
end
XM(t+1,:) = f(:,t);
PM(:,:,t+1) = A0(:,:,t)*PP(:,:,t)*A0(:,:,t)' + Q;
end
