function pegel = fun_kin_Welle_routing(pegel,pegel_no,vorgaenger,area)


% Pegelnummern der Pegel oberstrom (von denen jeweils geroutet wird):
p_vor_alle = vorgaenger;
% Anzahl der Pegelvorgänger (Laufvariable für Schleife):
p_no = 0;

% im Teilgebiet generierter Zufluss
Qges_short = pegel(1, pegel_no).discharge_generation(:,2)/area*3.6;

% Erweiterung in Minuten-Zeitschritte
Qges_MAT = repmat(Qges_short,1,60);
Qges_long = reshape(Qges_MAT',[],1)';

% Ausfluss des Routingabschnitts (hierzu werden später noch
% die von Oberliegern gerouteten Abflüsse addiert):
Q_out = Qges_long;
% aktueller Abfluss (der geroutet wird):
Qakt_long = nan(1,length(Qges_long));

for i = 1:length(p_vor_alle)
    p_vor = p_vor_alle(i);  % Routing von allen Vorgänger-Pegelpunkten
    p_no=p_no+1;
    
    % Parameter der Routing-Tabelle:
    intzahl = length(pegel(p_vor).routing.Q_in)-1;
    qmin = pegel(p_vor).routing.qmin;
    qmax = pegel(p_vor).routing.qmax;
    n = pegel(p_vor).routing.n;
    QvzuQ = pegel(p_vor).routing.QvzuQ;
    delta = pegel(p_vor).routing.delta;
    dq = log(qmax/qmin)/intzahl;
    subint = 60; % Minuten
    
    mi = subint + delta(1);
    fluss = nan(1,mi*n+subint);
    
    % Ganglinie des Vorgänger-Pegelpunkts (= Zuflussganglinie zum Flussabschnitt)
    if isfield(pegel(p_vor).routing, 'Qout') ==1
        Qakt_long_orig = pegel(p_vor).routing.Qout;
    else
        Qakt_short = pegel(p_vor).outflow(:,2)/pegel(p_vor).gage.area*3.6;
        % Erweiterung in Minuten-Zeitschritte
        Qakt_MAT = repmat(Qakt_short,1,60);
        Qakt_long_orig = reshape(Qakt_MAT',[],1)';
    end
    
    % Umrechnung von kh und kv
    kh = exp(-1/pegel(p_vor).routing.kh);
    kv = exp(-1/pegel(p_vor).routing.kv);
    
    t_short = 0;
    
    for t = 1:60:length(Qakt_long_orig)-59
        t_short = t_short+1;
        
        % aktueller Zeitschritt (60 min)
        Qakt = Qakt_long_orig(t:t+59);
        
        
        if t ==1
            Qaus = nan(n,length(pegel(p_vor).routing.Qaus(1,:)));
            qvzuq = nan(1,n);
            last_index = nan(1,n);
            
            for j = 1:n
                Qaus(j,:) = pegel(p_vor).routing.Qaus(j,:);
                qvzuq(j) = 0;
                last_index(j) = 360;
            end
            fluss = pegel(p_vor).routing.fluss';
        else
            fluss(last_index(n)+1:mi*n+subint) = Qakt(1);
        end
        

        for j = n:-1:1
            
            % Bereich in der W-Q-Tabelle
            %  [ber] = fun_getbereich(dq,max(Qakt),qmin,intzahl);
            q = max(Qakt)*277.7778;
            if q<0.0000001; q=0.0000001; end
            if dq>0
                ber = ceil(log(q/qmin)/dq);
            else
                ber=1;
            end
            if ber<1
                ber=1;
            elseif ber>intzahl+1
                ber=intzahl+1;
            end
            
            % Verschiebung in Minuten
            ver_zeitschritt = delta(ber);
            
            % Verschiebung auf Fluss-Array
            ver_fluss = subint + 1 + (j-1)*mi + ver_zeitschritt;
            
            if ver_fluss > last_index(j) +1
                
                % Lücken werden mit dem aktuellen (Qakt) oder
                % dem vorhergehenden (fluss) Abfluss gefüllt
                % (dem kleineren der Beiden)
                Qref = min(fluss(last_index(j)),Qakt(1));
                
                fluss(last_index(j)+1:ver_fluss-1) = Qref;
                
                % Kopieren des Qakt-Abflusses in den Fluss-array
                fluss(ver_fluss:ver_fluss+subint-1) = Qakt;%
            else
                % Kopieren des Qakt-Abflusses in den Fluss-array
                fluss(ver_fluss:ver_fluss+subint-1) = Qakt;
            end

            
            % Verschieben um konstant <subint>-Werte
            % (eigentliches Routing)
            fluss(1+(j-1)*mi:j*mi) = fluss(1+(j-1)*mi+subint:j*mi+subint);
            
            % neuen last_index berechnen
            last_index(j) = ver_fluss-1;
            
            % Bestimmen des Abflussbereichs des Ausflusses in
            % der W-Q-Tabelle
            fluss_max = max(fluss(1+(j-1)*subint:1+(j-1)*subint+subint-1));
            %  [ber] = fun_getbereich(dq,fluss_max,qmin,intzahl);
            q = fluss_max*277.7778;
            if q<0.0000001; q=0.0000001; end
            if dq>0
                ber = ceil(log(q/qmin)/dq);
            else
                ber=1;
            end
            if ber<1
                ber=1;
            elseif ber>intzahl+1
                ber=intzahl+1;
            end
            
            Q_v = fluss(1+(j-1)*mi:subint+(j-1)*mi).*QvzuQ(ber);
            Q_h = fluss(1+(j-1)*mi:subint+(j-1)*mi).*(1-QvzuQ(ber));
            
            Q_v_alt = Qaus(j,:).*qvzuq(j);
            Q_h_alt = Qaus(j,:).*(1-qvzuq(j));
            
            Qakt = Q_v_alt*kv + (Q_v*(1-kv)) + Q_h_alt*kh + (Q_h*(1-kh));
            t_range = t:t+59;
            Qakt_long(t_range) = Qakt;
            
            Qaus(j,:) = Qakt;
            qvzuq(j) = QvzuQ(ber);
        end
        
    end
    
    Q_out = Q_out + Qakt_long*pegel(pegel_no).Aezg(p_vor)/pegel(pegel_no).Aezg(pegel_no);
    
end
pegel(pegel_no).routing.Qout = Q_out;
Q_out_m3s = Q_out*pegel(pegel_no).gage.area/3.6;

Q_out_MAT = reshape(Q_out_m3s,60,length(Qakt_long_orig)/60)';
Q_out_ges = sum(Q_out_MAT,2)/60;

% Abspeichern der Daten in "pegel"
pegel(1, pegel_no).discharge_routed(:,2) = Q_out_ges - pegel(1, pegel_no).discharge_generation(:,2);                  % routed discharge
pegel(1, pegel_no).inflow(:,2) = Q_out_ges;

pegel(1, pegel_no).inflow(:,1) = pegel(1, pegel_no).inflow_wasim(:,1);  % date

end
