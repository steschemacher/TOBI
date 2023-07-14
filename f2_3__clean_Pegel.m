function [pegel] = f2_3__clean_Pegel(pegel,E_no,Ereignis,richards,max_gage)

% Abflüsse einfügen
for i = 1:max_gage
    pegel(i).pegel_addAbfluss(i,Ereignis(E_no).qgko_abs,Ereignis(E_no).qgko_date,[],Ereignis(E_no).gwin_abs,richards(E_no,:),...
        Ereignis(E_no).qdir_abs,Ereignis(E_no).qifl_abs,Ereignis(E_no).qbas_abs);
end

% Ganglinien aus vorherigen Simulationen löschen
for i = 1:max_gage
    pegel(i).inflow = [];
    pegel(i).outflow = [];
    pegel(i).inflow{1} = pegel(i).inflow_wasim;      % no upstream gage -> no routing required
    pegel(i).outflow = pegel(i).inflow;           % no basin at the gage -> no basin retention
    pegel(i).discharge_routed(:,1) = pegel(i).inflow{1}(:,1);    % time
%     pegel(i).discharge_routed(:,2) = pegel(i).inflow{1}(:,2)-pegel(i).discharge_generation(:,2);  % routed discharge
    
    pegel(i).discharge_routed(:,2) = pegel(i).inflow{1}(:,2)-pegel(i).discharge_qdir(:,2)...
        -pegel(i).discharge_qifl(:,2)-pegel(i).discharge_qbas(:,2);  % routed discharge
    pegel(i).neighbors.becken_vorgaenger = [];
    pegel(i).becken_vor_yn = 0;
    pegel(i).becken_yn = 0;
    pegel(i).beckenparameter = [];
    pegel(i).becken = [];
%     pegel(i).abminderung = [];
    pegel(i).VQ_fest_yn = 0;
    pegel(i).abstraction_yn = 0;
    pegel(i).gage.S_gesamt = [0,0,0];
end