function pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,k)
            %---------------------------------------------------------------%
            % function to calculate the basin retention (determine outflow) %
            %---------------------------------------------------------------%
%             if exist('k','var')==0
%                 k=1;
%             end
%             if isfield(pegel(pegel_no).beckenparameter(k), 'difference')==1
%                 difference = pegel(pegel_no).beckenparameter(k).difference;
%             end
            
            dt = 3600;       % time step for basin retention (time step of the input data)
            
            
            % definition of the time stamps of the variables
%             pegel(pegel_no).becken.h(:,1) = pegel(pegel_no).inflow(:,1);

%             pegel(pegel_no).becken(k).S(:,1) = inflow_tmp(:,1);
%             pegel(pegel_no).outflow(:,1) = inflow_tmp(:,1);
            
%             pegel(pegel_no).becken(k).S(1,2) = pegel(pegel_no).beckenparameter(k).S0;               % 1st step: basin is empty
%             pegel(pegel_no).outflow(1,2) = pegel(pegel_no).inflow(1,2);  % 1st step: outflow = inflow (assumption that inflow < baseflow)
            if k==1
                inflow_tmp = pegel(pegel_no).inflow{k};
            else
                inflow_tmp = pegel(pegel_no).outflow{k-1};
            end
            
            laenge = length(inflow_tmp(:,2));    % time series
            
            outflow_tmp(1:length(inflow_tmp(:,1)),1) = inflow_tmp(:,1);
            outflow_tmp(1,2) = inflow_tmp(1,2);
            S(1:length(inflow_tmp(:,1)),1) = inflow_tmp(:,1);
            S(1,2) = pegel(pegel_no).beckenparameter(k).S0;
            Q_maxF = pegel(pegel_no).beckenparameter(k).Q_maxF;
%             S = pegel(pegel_no).becken(k).S;
            ShQ = pegel(pegel_no).beckenparameter(k).ShQ;
            for t= 2:laenge
                % if the basin is empty and the inflow can be discharged in free flow:
                if (S(t-1,2)<=0 && inflow_tmp(t,2)<Q_maxF)==1
                    outflow_tmp(t,2) = inflow_tmp(t,2);
                % basin retention (storage equation):
                else
                    % known values of the storage equation -> left side
                    % (inflow, outflow + storage at the last time step):
                    bekannt = .5*dt*(inflow_tmp(t,2)+inflow_tmp(t-1,2)-outflow_tmp(t-1,2))+S(t-1,2);
                    % lower and upper limit of the V-Q-relationship:
                    I_links = 1;
                    I_rechts = length(ShQ);
                    
%                     begrenzung = [I_links,I_rechts];          % limits
                    I_mitte = round((I_links+I_rechts)/2);    % value in the middle
                    % loop while lower and upper limits are no direct neighbors:
                    while (I_links~=I_mitte && I_rechts~=I_mitte)==1%any(begrenzung==I_mitte)==0%ismember(begrenzung,I_mitte)==0
                        % unknown values of the storage equation -> right side
                        % -> values of the middle between both limits:
                        gesucht_mitte = ShQ(I_mitte,1)+.5*dt*ShQ(I_mitte,3);
                        % procedure to find next values for the iteration:
                        if gesucht_mitte > bekannt    % middle is too large
                            I_rechts = I_mitte;       % right limit = middle
                        else                          % middle is too small
                            I_links = I_mitte;        % left limit  = middle
                        end
                        I_mitte = round((I_links+I_rechts)/2);    % find new middle
%                         begrenzung = [I_links,I_rechts];          % define new limits
                    end
                    % result: storage-discharge combination lies in between
                    % two points of the relation -> linear interpolation
                    % gradient of the line between both points
                    m = (ShQ(I_rechts,1)-ShQ(I_links,1))/(ShQ(I_rechts,3)-ShQ(I_links,3));
                    % distance between left discharge and wanted discharge:
                    x = (bekannt - ShQ(I_links,1) - .5*dt*ShQ(I_links,3))/(m+.5*dt);
                    % interpolated discharge and storage:
                    outflow_tmp(t,2) = ShQ(I_links,3)+x;
                    S(t,2) = ShQ(I_links,1)+m*x;
                   
%                     if exist('difference')==1
%                         outflow_tmp(t,2) = outflow_tmp(t,2) - difference(t);
%                         S(t,2) = S(t,2) + difference(t)*3600;
%                     end
                    
                    
                    % conditions when inflow = outflow and S = 0:
                    % (inflow hydrograph falling + outflow hydrograph falling + inflow > outflow + S<100*outflow)
                    if outflow_tmp(t-1,2)>outflow_tmp(t,2) && inflow_tmp(t-1,2)>inflow_tmp(t,2) && outflow_tmp(t,2)<inflow_tmp(t,2) && S(t,2)<50000
                        outflow_tmp(t,2) = inflow_tmp(t,2);
                        S(t,2)=0;
                    end
                    if outflow_tmp(t,2)<0
                        outflow_tmp(t,2) = inflow_tmp(t,2);
                        S(t,2)=0;
                    end
                    
                    if isfield(pegel(pegel_no).beckenparameter(k),'S_max')
                    dt_neu = dt;
                    if S(t,2) < 1.1*outflow_tmp(t,2)*dt_neu && S(t,2) > pegel(pegel_no).beckenparameter(k).S_max
                    
                    while S(t,2) < 1.1*outflow_tmp(t,2)*dt_neu && S(t,2) > pegel(pegel_no).beckenparameter(k).S_max && dt_neu>200
%                         teiler = ceil(round(S(t,2)/outflow_tmp(t,2)))
                        dt_neu = dt_neu/2;
                        
                        in_tmp_vor = inflow_tmp(t-1,2);
                        out_tmp_vor = outflow_tmp(t-1,2);
                        speicher_tmp_vor = S(t-1,2);
                        
                        for schritt = 1:dt/dt_neu
                            in_tmp = (inflow_tmp(t-1,2)*(dt/dt_neu-schritt)+inflow_tmp(t,2)*schritt)/(dt/dt_neu);
                            
                            % known values of the storage equation -> left side
                            % (inflow, outflow + storage at the last time step):
                            bekannt = .5*dt_neu*(in_tmp+in_tmp_vor-out_tmp_vor)+speicher_tmp_vor;
                            % lower and upper limit of the V-Q-relationship:
                            I_links = 1;
                            I_rechts = length(ShQ);
                            
                            %                     begrenzung = [I_links,I_rechts];          % limits
                            I_mitte = round((I_links+I_rechts)/2);    % value in the middle
                            % loop while lower and upper limits are no direct neighbors:
                            while (I_links~=I_mitte && I_rechts~=I_mitte)==1%any(begrenzung==I_mitte)==0%ismember(begrenzung,I_mitte)==0
                                % unknown values of the storage equation -> right side
                                % -> values of the middle between both limits:
                                gesucht_mitte = ShQ(I_mitte,1)+.5*dt_neu*ShQ(I_mitte,3);
                                % procedure to find next values for the iteration:
                                if gesucht_mitte > bekannt    % middle is too large
                                    I_rechts = I_mitte;       % right limit = middle
                                else                          % middle is too small
                                    I_links = I_mitte;        % left limit  = middle
                                end
                                I_mitte = round((I_links+I_rechts)/2);    % find new middle
                                %                         begrenzung = [I_links,I_rechts];          % define new limits
                            end
                            % result: storage-discharge combination lies in between
                            % two points of the relation -> linear interpolation
                            
                            % gradient of the line between both points
                            m = (ShQ(I_rechts,1)-ShQ(I_links,1))/(ShQ(I_rechts,3)-ShQ(I_links,3));
                            % distance between left discharge and wanted discharge:
                            x = (bekannt - ShQ(I_links,1) - .5*dt_neu*ShQ(I_links,3))/(m+.5*dt_neu);
                            % interpolated discharge and storage:
                            out_tmp = ShQ(I_links,3)+x;
                            speicher_tmp = ShQ(I_links,1)+m*x;
                            
                            in_tmp_vor = in_tmp;
                            out_tmp_vor = out_tmp;
                            speicher_tmp_vor = speicher_tmp;
  
                        end
                        
                        
                    end
                    
                    outflow_tmp(t,2) = out_tmp;
                    S(t,2) = speicher_tmp;
                    
%                     if exist('difference')==1
%                         outflow_tmp(t,2) = outflow_tmp(t,2) - difference(t);
%                         S(t,2) = S(t,2) + difference(t)*3600;
%                     end
                    
                    % conditions when inflow = outflow and S = 0:
                    % (inflow hydrograph falling + outflow hydrograph falling + inflow > outflow + S<100*outflow)
%                     if pegel(pegel_no).gage.no==45
%                     end
                    if outflow_tmp(t-1,2)>outflow_tmp(t,2) && inflow_tmp(t-1,2)>inflow_tmp(t,2) && outflow_tmp(t,2)<inflow_tmp(t,2) && S(t,2)<50000
                        outflow_tmp(t,2) = inflow_tmp(t,2);
                        S(t,2)=0;
                    end
                    if outflow_tmp(t,2)<0
                        outflow_tmp(t,2) = inflow_tmp(t,2);
                        S(t,2)=0;
                    end
                    end
                    end
%                     pegel(pegel_no).becken(k).S_max = max(pegel(pegel_no).becken(k).S(:,2));
                end
%                 pegel(pegel_no).becken(k).S_max = max(pegel(pegel_no).becken(k).S(:,2));
            end
            pegel(pegel_no).outflow{k} = outflow_tmp;
            pegel(pegel_no).inflow{k} = inflow_tmp;
            pegel(pegel_no).becken(k).S = S;
            pegel(pegel_no).becken(k).S_max = max(pegel(pegel_no).becken(k).S(:,2));
            if isfield(pegel(pegel_no).routing, 'Qout') ==1
                pegel(pegel_no).routing = rmfield(pegel(pegel_no).routing,'Qout');
            end
            
        end