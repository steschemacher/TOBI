function f2_5__retention_routing(pegel,pegel_no,method)
% pegel_no
%---------------------------------------------------------------%
%  function to route the discharges within the basin            %
%  (recursive fuction which starts at one point and calculates  % 
%    all missing hydrographs)                                   %
%---------------------------------------------------------------%

% if the end of the routing is reached 
% (last gage didn't have a downstream neighbor)
if isempty(pegel_no)==1
    return
end
   
%% gage with basin and without upstream basin (= potential start gage)
if pegel(pegel_no).becken_yn ==1 && pegel(pegel_no).becken_vor_yn ==0
    pegel(pegel_no).inflow{1} = pegel(pegel_no).inflow_wasim;      % inflow = simulated discharge (WaSiM)  
    if pegel(pegel_no).VQ_fest_yn == 0     % V-h über 
        f2_5_1__retention_optimierung(pegel,pegel_no,method);

    else
        pegel(pegel_no).inflow{1}(:,2) = pegel(pegel_no).inflow{1}(:,2)+pegel(pegel_no).beckenparameter.difference(:,1);
%         pegel(pegel_no).retention;
        pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,1);

    end
    
    % abstraction at the gage
    if pegel(pegel_no).abstraction_yn == 1
        pegel(pegel_no).abstraction;
    end
    
    pegel = fun_S_gesamt(pegel,pegel_no);

    f2_5__retention_routing(pegel,pegel(pegel_no).neighbors.nachfolger,method);    % restart function with next gage downstream 
    return
end

%% gage with upstream basin (no potential start gage)
%%% gage with ONE upstream neighbor
if length(pegel(pegel_no).neighbors.vorgaenger)<2
    % no additional inflows need to be regarded -> only routing from one gage to the other)
    
    if isempty(pegel(1, pegel(pegel_no).neighbors.vorgaenger(1)).outflow) ==0
        % if the infow of the upstream neighbor has already been defined
        
%         pegel = muskingum(pegel,pegel_no,[]);              % inflow = routed outflow from upstream neighbor
%         pegel = fun_kin_Welle_routing(pegel,pegel_nopegel(pegel_no).neighbors.vorgaenger,...
%             pegel(pegel_no).gage.area);              % inflow = routed outflow from upstream neighbor
vorgaenger = pegel(1, pegel_no).neighbors.vorgaenger;
% discharge_generation = pegel(pegel_no).discharge_generation;
discharge_generation = pegel(pegel_no).discharge_qdir+pegel(pegel_no).discharge_qifl+pegel(pegel_no).discharge_qbas;
pno_area = pegel(pegel_no).gage.area;
Aezg = pegel(pegel_no).Aezg;
discharge_share_comb = zeros(size(discharge_generation,1),3);
for p_vor = 1:length(vorgaenger)    
    area(p_vor) = pegel(vorgaenger(p_vor)).gage.area;
    intzahl(p_vor) = length(pegel(vorgaenger(p_vor)).routing.Q_in)-1;
    qmin(p_vor) = pegel(vorgaenger(p_vor)).routing.qmin;
    qmax(p_vor) = pegel(vorgaenger(p_vor)).routing.qmax;
    n(p_vor) = pegel(vorgaenger(p_vor)).routing.n;
    QvzuQ(p_vor,:) = pegel(vorgaenger(p_vor)).routing.QvzuQ;
    delta(p_vor,:) = pegel(vorgaenger(p_vor)).routing.delta;
    if isfield(pegel(vorgaenger(p_vor)).routing, 'Qout') ==1
        Qout{p_vor} = pegel(vorgaenger(p_vor)).routing.Qout;
    else
        Qout{p_vor} = nan(2,2);
    end
    outflow(p_vor,:) = pegel(vorgaenger(p_vor)).outflow{end}(:,2)';
    routing_kh(p_vor) = pegel(vorgaenger(p_vor)).routing.kh;
    routing_kv(p_vor) = pegel(vorgaenger(p_vor)).routing.kv;
    Qaus_input{p_vor} = pegel(vorgaenger(p_vor)).routing.Qaus;
    routing_fluss{p_vor} = pegel(vorgaenger(p_vor)).routing.fluss;
    discharge_share{p_vor} = pegel(vorgaenger(p_vor)).discharge_shares*area(p_vor);
    discharge_share_comb = discharge_share_comb + discharge_share{p_vor};
end
[pegel(pegel_no).routing.Qout,pegel(1, pegel_no).discharge_routed(:,2),...
    pegel(1, pegel_no).inflow] = f2_5_3__routing_kin_welle...
    (pegel_no,vorgaenger,discharge_generation,pno_area,area,intzahl,qmin,qmax,...
    n,QvzuQ,delta,Qout,outflow,routing_kh,routing_kv,Qaus_input,routing_fluss,...
    Aezg);
discharge_total = pegel(pegel_no).discharge_qdir(:,2)+pegel(pegel_no).discharge_qifl(:,2)+pegel(pegel_no).discharge_qbas(:,2);
pegel_discharge_shares = [pegel(pegel_no).discharge_qdir(:,2)./discharge_total,...
    pegel(pegel_no).discharge_qifl(:,2)./discharge_total,pegel(pegel_no).discharge_qbas(:,2)./discharge_total];
pegel(pegel_no).discharge_shares = (pegel_discharge_shares.*(pno_area-sum(area)) + discharge_share_comb)./pno_area;



        if pegel(pegel_no).becken_yn ==1
            % if there is a retention basin at the gage
            if pegel(pegel_no).VQ_fest_yn == 0
                f2_5_1__retention_optimierung(pegel,pegel_no,method);

%                 pegel(regional).becken.S_gesamt = S_gesamt+pegel(pegel_no).becken.S_max;
%             pegel(pegel_no).volumen_abfluss_beziehung;  % definition of the V-Q-relation of the basin
%             pegel(pegel_no).retention;                  % outflow = discharge after basin retention
            else
                pegel(pegel_no).inflow{1}(:,2) = pegel(pegel_no).inflow{1}(:,2)+pegel(pegel_no).beckenparameter.difference(:,1);
%                 pegel(pegel_no).retention;
                pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,1);

            end           
        else
             % if not
            pegel(pegel_no).outflow = pegel(pegel_no).inflow;   % outflow = inflow (no changes at the gage)
        end
        % abstraction at the gage
        if pegel(pegel_no).abstraction_yn == 1
            pegel(pegel_no).abstraction;
        end
        % find downstream neighbor and restart the function:
        pegel = fun_S_gesamt(pegel,pegel_no);
 
        f2_5__retention_routing(pegel,pegel(pegel_no).neighbors.nachfolger,method);
        return  
    else
        % if inflow of upstream neighbor needs to be calculated -> restart the function with that gage number: 
        if pegel(1, pegel(pegel_no).neighbors.vorgaenger(1)).becken_vor_yn ==1
            % if the upstream neighbor has a basin upstream -> restart the function with the number of that gage
            pegel = fun_S_gesamt(pegel,pegel_no);
 
            f2_5__retention_routing(pegel,pegel(pegel(pegel_no).neighbors.vorgaenger(1)).neighbors.becken_vorgaenger(1),method);
            return
        elseif pegel(1, pegel(pegel_no).neighbors.vorgaenger(1)).becken_yn ==1
            % if the upstream neighbor has NO basin upstream, but there is a basin located at the gage itself
            pegel = fun_S_gesamt(pegel,pegel_no);
 
            f2_5__retention_routing(pegel,pegel(pegel_no).neighbors.vorgaenger(1),method);
            return
        end
    end
    
%%% gage with MORE THAN ONE upstream neighbor (after inflow)
else
    % loop to check all upstream neighbors:
    for i = 1:length(pegel(pegel_no).neighbors.vorgaenger)        
        if isempty(pegel(1, pegel(pegel_no).neighbors.vorgaenger(i)).outflow) ==1
            % if there is NO inflow calculated at the upstream neighbor       
            if pegel(1, pegel(pegel_no).neighbors.vorgaenger(i)).becken_vor_yn ==1
                % if the upstream neighbor has a basin upstream -> restart the function with the number of that gage
                pegel = fun_S_gesamt(pegel,pegel_no);
 
                f2_5__retention_routing(pegel,pegel(pegel(pegel_no).neighbors.vorgaenger(i)).neighbors.becken_vorgaenger(1),method);       
                return
            elseif pegel(1, pegel(pegel_no).neighbors.vorgaenger(i)).becken_yn ==1
                % if the upstream neighbor has NO basin upstream, but there is a basin located at the gage itself
                pegel = fun_S_gesamt(pegel,pegel_no);
 
                f2_5__retention_routing(pegel,pegel(pegel_no).neighbors.vorgaenger(i),method);
                return
            end
            return
        end
    end
    

    % superposition of upstream hydrographs to calculate new hydrograph of the
    % gage (including routing)
% 
        QA_geg = 0;
        Qzu_geg = 0;
        for i = 1:length(pegel(pegel_no).neighbors.vorgaenger)
            QA_geg = QA_geg +  pegel(pegel(pegel_no).neighbors.vorgaenger(i)).inflow_wasim(:,3);  % original hydrograph (neighbor)
            Qzu_geg = Qzu_geg + pegel(pegel(pegel_no).neighbors.vorgaenger(i)).outflow{end}(:,2);          % new routed hydrograph (neighbor)
        end
        zufluesse = struct('QA_geg',QA_geg,'Qzu_geg',Qzu_geg);
        pegel(pegel_no).becken.zufluesse = zufluesse;
%         pegel = muskingum(pegel,pegel_no,zufluesse);              % inflow = routed outflow from upstream neighbor
%         pegel = fun_kin_Welle_routing(pegel,pegel_no,pegel(pegel_no).neighbors.vorgaenger,...
%             pegel(pegel_no).gage.area);              % inflow = routed outflow from upstream neighbor
vorgaenger = pegel(1, pegel_no).neighbors.vorgaenger;
% discharge_generation = pegel(pegel_no).discharge_generation;
discharge_generation = pegel(pegel_no).discharge_qdir+pegel(pegel_no).discharge_qifl+pegel(pegel_no).discharge_qbas;
pno_area = pegel(pegel_no).gage.area;
Aezg = pegel(pegel_no).Aezg;
discharge_share_comb = zeros(size(discharge_generation,1),3);
for p_vor = 1:length(vorgaenger)
    area(p_vor) = pegel(vorgaenger(p_vor)).gage.area;
    intzahl(p_vor) = length(pegel(vorgaenger(p_vor)).routing.Q_in)-1;
    qmin(p_vor) = pegel(vorgaenger(p_vor)).routing.qmin;
    qmax(p_vor) = pegel(vorgaenger(p_vor)).routing.qmax;
    n(p_vor) = pegel(vorgaenger(p_vor)).routing.n;
    QvzuQ(p_vor,:) = pegel(vorgaenger(p_vor)).routing.QvzuQ;
    delta(p_vor,:) = pegel(vorgaenger(p_vor)).routing.delta;
    if isfield(pegel(vorgaenger(p_vor)).routing, 'Qout') ==1
        Qout{p_vor} = pegel(vorgaenger(p_vor)).routing.Qout;
    else
        Qout{p_vor} = nan(2,2);
    end
    outflow(p_vor,:) = pegel(vorgaenger(p_vor)).outflow{end}(:,2)';
    routing_kh(p_vor) = pegel(vorgaenger(p_vor)).routing.kh;
    routing_kv(p_vor) = pegel(vorgaenger(p_vor)).routing.kv;
    Qaus_input{p_vor} = pegel(vorgaenger(p_vor)).routing.Qaus;
    routing_fluss{p_vor} = pegel(vorgaenger(p_vor)).routing.fluss;
    discharge_share{p_vor} = pegel(vorgaenger(p_vor)).discharge_shares*area(p_vor);
    discharge_share_comb = discharge_share_comb + discharge_share{p_vor};
end
[pegel(pegel_no).routing.Qout,pegel(1, pegel_no).discharge_routed(:,2),...
    pegel(1, pegel_no).inflow] = f2_5_3__routing_kin_welle...
    (pegel_no,vorgaenger,discharge_generation,pno_area,area,intzahl,qmin,qmax,...
    n,QvzuQ,delta,Qout,outflow,routing_kh,routing_kv,Qaus_input,routing_fluss,...
    Aezg);
discharge_total = pegel(pegel_no).discharge_qdir(:,2)+pegel(pegel_no).discharge_qifl(:,2)+pegel(pegel_no).discharge_qbas(:,2);
pegel_discharge_shares = [pegel(pegel_no).discharge_qdir(:,2)./discharge_total,...
    pegel(pegel_no).discharge_qifl(:,2)./discharge_total,pegel(pegel_no).discharge_qbas(:,2)./discharge_total];
pegel(pegel_no).discharge_shares = (pegel_discharge_shares.*(pno_area-sum(area)) + discharge_share_comb)./pno_area;


  
    % routed discharge
%     pegel(pegel_no).discharge_routed(:,2) = pegel(pegel_no).inflow{1}(:,2)-pegel(pegel_no).discharge_generation(:,2);
    pegel(pegel_no).discharge_routed(:,2) = pegel(pegel_no).inflow{1}(:,2)...
        -pegel(pegel_no).discharge_qdir(:,2)-pegel(pegel_no).discharge_qifl(:,2)-pegel(pegel_no).discharge_qbas(:,2);

    
    if pegel(pegel_no).becken_yn ==0    % no basin at the gage
        pegel(pegel_no).outflow = pegel(pegel_no).inflow;
    else                                % basin at the gage
        if pegel(pegel_no).VQ_fest_yn == 0
            % V-Q-relation + basin retention
           f2_5_1__retention_optimierung(pegel,pegel_no,method);

        else
            pegel(pegel_no).inflow{1}(:,2) = pegel(pegel_no).inflow{1}(:,2)+pegel(pegel_no).beckenparameter.difference(:,1);
%             pegel(pegel_no).retention;
            pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,1);
            
        end
    end
    % abstraction at the gage
    if pegel(pegel_no).abstraction_yn == 1
        pegel(pegel_no).abstraction;
    end
    % restart function with downstream neighbor
    pegel = fun_S_gesamt(pegel,pegel_no);
 
    f2_5__retention_routing(pegel,pegel(pegel_no).neighbors.nachfolger,method);

return
end

end