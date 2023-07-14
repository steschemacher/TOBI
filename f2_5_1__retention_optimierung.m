function f2_5_1__retention_optimierung(pegel,pegel_no,method)
%------------------------------------------------------------%
%  function to determine the basin retention                 %
%  (including optimization of the basin opening dimension)  %
%------------------------------------------------------------%


if pegel(pegel_no).beckenparameter(1).S_max==0
    pegel(pegel_no).outflow{1} = pegel(pegel_no).inflow{1};
else
    hoehen = zeros(1,2);
    for i = 1:length(pegel(pegel_no).beckenparameter)
        hoehen(i,2) = pegel(pegel_no).beckenparameter(i).geo;
    end
    
    % if length(hoehen(:,1))>1
    hoehen = sortrows(hoehen,-2);
    hoehen(:,1) = 1:length(hoehen(:,2));
    % end
    AL_rechts = zeros(1,1);
    for k = 1:max(hoehen(:,1))
        if method ~= 3
            % set upper and lower limits for the opening dimension
            if length(pegel(pegel_no).gage.AL_min)>1
                AL_links = pegel(pegel_no).gage.AL_min(k);   % lower limit
                AL_rechts = min(pegel(pegel_no).gage.AL_max(k),pegel(pegel_no).beckenparameter(k).h_W^2) ; %(max(pegel(pegel_no).inflow(:,2))^2/(pegel(pegel_no).beckenparameter.mu_A^2*9.81))^(2/5); % upper limit
            else
                AL_links = pegel(pegel_no).gage.AL_min;   % lower limit
                AL_rechts = min(pegel(pegel_no).gage.AL_max,pegel(pegel_no).beckenparameter(k).h_W^2) ; %(max(pegel(pegel_no).inflow(:,2))^2/(pegel(pegel_no).beckenparameter.mu_A^2*9.81))^(2/5); % upper limit
            end
            
            %     if k>1
            %         pegel(pegel_no).inflow{k} = pegel(pegel_no).outflow{k};
            %     end
            % check boundary conditions
            % storage exceeds maximum storage at AL_rechts
            pegel(pegel_no).beckenparameter(k).A_L = AL_rechts;
            if isfield(pegel(pegel_no).beckenparameter(k), 'ShAs')==1
                volumen_abfluss_beziehung_real(pegel(pegel_no),method,k);	% generation of the V-Q-relation of the retention basin
            else
                volumen_abfluss_beziehung(pegel(pegel_no),method);	% generation of the V-Q-relation of the retention basin
            end
            % pegel(pegel_no).retention(k);               	% outflow = discharge after basin retention calculation
            pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,k);
            maximales_S = max(pegel(pegel_no).becken(k).S(:,2));   % maximum fill of the storage for chosen A_L
            if maximales_S > pegel(pegel_no).beckenparameter(k).S_max
                %     fprintf('Hochwasserentlastung von Becken %d ist notwendig. Maximales Speichervolumen: %4.1f m^3\n',pegel_no,maximales_S)
                %     fprintf('Öffnungsgröße: %4.2f m^2\n',pegel(pegel_no).beckenparameter.A_L)
                pegel(pegel_no).beckenfuellung(k) = maximales_S/pegel(pegel_no).beckenparameter(k).S_max;
                AL_links = [];
                AL_rechts = [];
                maximales_S = [];
                %     clear AL_links AL_rechts maximales_S
                continue
                %     return
            end
            
            % delete values of previous runs
            %     if k==1
            pegel(pegel_no).beckenparameter(k).ShQ = [];
            pegel(pegel_no).becken(k) = [];
            pegel(pegel_no).outflow{k} = [];
            %     end
            
            % basin is not filled at AL_links
            pegel(pegel_no).beckenparameter(k).A_L = AL_links;
            if isfield(pegel(pegel_no).beckenparameter(k), 'ShAs')==1
                volumen_abfluss_beziehung_real(pegel(pegel_no),method,k);	% generation of the V-Q-relation of the retention basin
            else
                volumen_abfluss_beziehung(pegel(pegel_no),method);	% generation of the V-Q-relation of the retention basin
            end
            % pegel(pegel_no).retention(k);               	% outflow = discharge after basin retention calculation
            pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,k);
            maximales_S = max(pegel(pegel_no).becken(k).S(:,2));   % maximum fill of the storage for chosen A_L
            if maximales_S < pegel(pegel_no).beckenparameter(k).S_max;
                %     fprintf('Becken %d wird nicht gefüllt. Maximales Speichervolumen: %4.1f m^3\n',pegel_no,maximales_S)
                %     fprintf('Öffnungsgröße: %4.2f m^2\n',pegel(pegel_no).beckenparameter.A_L)
                pegel(pegel_no).beckenfuellung(k) = maximales_S/pegel(pegel_no).beckenparameter(k).S_max;
                AL_links = [];
                AL_rechts = [];
                maximales_S = [];
                %     clear AL_links AL_rechts maximales_S
                continue
                %     return
            end
            
            % delete values of previous runs
            pegel(pegel_no).beckenparameter(k).ShQ = [];
            pegel(pegel_no).becken(k) = [];
            pegel(pegel_no).outflow{k} = [];
            
            
            % start conditions: basin retention with mean opening size
            pegel(pegel_no).beckenparameter(k).A_L = (AL_rechts+AL_links)/2;
            if isfield(pegel(pegel_no).beckenparameter(k), 'ShAs')==1
                volumen_abfluss_beziehung_real(pegel(pegel_no),method,k);	% generation of the V-Q-relation of the retention basin
            else
                volumen_abfluss_beziehung(pegel(pegel_no),method);	% generation of the V-Q-relation of the retention basin
            end
            % pegel(pegel_no).retention(k);               	% outflow = discharge after basin retention calculation
            pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,k);
            maximales_S = max(pegel(pegel_no).becken(k).S(:,2));   % maximum fill of the storage for chosen A_L
            
            maximales_S_alt = 100000;
            % loop until maximum fill equals the defined storage volume
            while abs(maximales_S-pegel(pegel_no).beckenparameter(k).S_max)>1 && abs(maximales_S_alt-maximales_S)>1
                if maximales_S~=0
                    maximales_S_alt = maximales_S;
                else
                    maximales_S_alt = 100000;
                end
                
                if maximales_S < pegel(pegel_no).beckenparameter(k).S_max % basin is not full (opening should be smaller)
                    AL_rechts = pegel(pegel_no).beckenparameter(k).A_L;
                else    % emergency spillway is used (opening should be larger)
                    AL_links = pegel(pegel_no).beckenparameter(k).A_L;
                end
                
                % delete values of previous runs
                pegel(pegel_no).beckenparameter(k).ShQ = [];
                pegel(pegel_no).becken(k) = [];
                pegel(pegel_no).outflow{k} = [];
                
                pegel(pegel_no).beckenparameter(k).A_L = (AL_rechts+AL_links)/2;
                
                if isfield(pegel(pegel_no).beckenparameter(k), 'ShAs')==1
                    volumen_abfluss_beziehung_real(pegel(pegel_no),method,k);	% generation of the V-Q-relation of the retention basin
                else
                    volumen_abfluss_beziehung(pegel(pegel_no),method);	% generation of the V-Q-relation of the retention basin
                end
                %     pegel(pegel_no).retention(k);                 % outflow = discharge after basin retention calculation
                pegel = f2_5_2__retention_speichergleichung(pegel,pegel_no,k);
                maximales_S = max(pegel(pegel_no).becken(k).S(:,2));   % maximum fill of the storage for chosen A_L
                
            end
            % fprintf('Becken %d wird gefüllt. Maximales Speichervolumen: %4.1f m^3\n',pegel_no,maximales_S)
            % fprintf('Öffnungsgröße: %4.2f m^2\n',pegel(pegel_no).beckenparameter.A_L)
            pegel(pegel_no).beckenfuellung(k) = maximales_S/pegel(pegel_no).beckenparameter(k).S_max;
            AL_links = [];
            AL_rechts = [];
            maximales_S = [];
            % clear AL_links AL_rechts maximales_S
            
        else
            pegel = f2_5_4__retention_gesteuert(pegel,pegel_no,k);
        end
    end
end

return