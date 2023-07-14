classdef Pegel<handle
    
    properties
        % general variables
        river           % struct: river no, flowpoints, gages
        riverpoints     % table: flowpoint no, col, row, river no, gage no
        grids           % DGM, RIV, DEP, EZG, FZS, FLK, PUR
        Aezg            % size of the subcatchments [km^2]
        Aezg_share
        
        % hydrographs
        inflow_wasim            % total discharge: qgko (WaSiM)
        discharge_generation    % generated discharge: qges (WaSiM)
        discharge_qdir
        discharge_qifl
        discharge_qbas
        discharge_shares
        discharge_routed        % routed discharge (Matlab: Muskingum)
        gw_infiltration         % infiltration from rivers to the groundwater: gwin (WaSiM)
        inflow                  % inflow to the gage (after routing, before basin)
        outflow                 % outflow from the gage (after basin)
        routing
        %         fluss
        %         qvzuq
        
        % location, neighbors and properties
        neighbors           % struct: vorgaenger, nachfolger, alle_nachfolger, becken_vorgaenger
        gage                % struct: no, river, area, coordinates (col, row), elevation
        beckenparameter     % struct: HQ100,BHQ,Q_maxF,Q_maxL,Q_maxG, S_max,h_W,b_W,h_D,f,m,a, mu_A,mu_W,A_L,AL_min, ShQ (table)
        becken              % struct: h, S
        abstr_par           % struct: min_abstr, max_abtr, frac_abstr
        pegel_vor_yn        % gage upstream of that gate: 1/0 (yes/no)
        
        % variables for gage selection
        becken_moeglich     % possibility for basin: 1/0 (yes/no)
        becken_yn           % basin at that gage: 1/0 (yes/no)
        becken_vor_yn       % basin upstream of that gage: 1/0 (yes/no)
        VQ_fest_yn          % volume-discharge-relation: 1/0 (yes/no)
        abstraction_yn      % abstraction rule: 1/0 (yes/no)
        beckenfuellung      % relative filling of the basin (compared to S_max)
        
        % result
        %         abminderung
                valid
    end
    
    methods
        function obj=Pegel(gage_no,grids,river,A_ezg,river_points,MQ,HHQ)
            %---------------------------------------------------------------%
            % function to set initial values to the variables of every gage %
            %---------------------------------------------------------------%
            
            %%% flow points, gages and neighbors
            % information about the whole river net or all gages
            obj.grids = grids;                  % DGM, RIV, DEP, EZG, FZS, FLK, PUR
            obj.river = river;                  % all rivers, flow depth, flow points and gages in the catchment
            obj.riverpoints = river_points;     % table with information about river points: [ID,x-coord,y-coord,river_ID,gage_ID]
            obj.Aezg = A_ezg;                   % areas of the subcatchments [km^2]
            obj.Aezg_share = A_ezg./A_ezg(1);
            
            % information about the gage
            obj.gage.no = gage_no;                                          % gage ID
            obj.gage.river = river_points(river_points(:,5)==gage_no,4);    % river ID
            obj.gage.area = A_ezg(1,gage_no);                               % subcatchment area [km^2]
            row = obj.river(obj.gage.river).gages(:,1)==obj.gage.no;
            obj.gage.coordinates = obj.river(obj.gage.river).gages(row,2:3);% coordinates [row,col]
            obj.gage.elevation = obj.river(obj.gage.river).gages(row,4);    % elevation [m+NN]
            
            obj.gage.S_gesamt = [0,0,0];
            
            % definition of initial values of variables
            obj.neighbors.nachfolger = [];          % next gage downstream
            obj.neighbors.vorgaenger = [];          % next gage upstream
            obj.neighbors.alle_nachfolger = [];     % all gages downstream
            obj.neighbors.becken_vorgaenger = [];   % all gages with basins upstream
            obj.becken_vor_yn = 0;                  % there is a basin upstream (yes/no = 1/0)
            obj.becken_yn = 0;                      % this gage has a basin (yes/no = 1/0)
            obj.routing = [];
            obj.pegel_vor_yn = 0;                   % there is a gage upstream (yes/no = 1/0)
            %             obj.becken_moeglich = 0;                % alle Beckenstandorte werden erstmal ausgeschlossen (später einfügen der real möglichen Standorte)
            
            %%% parameters for the basin design (design flood and discharge for dimensioning
            %             tageswert = reshape(qgko_abs(1:end-1,gage_no),[],24);
            obj.gage.MQ = MQ(obj.gage.no);%MQ/obj.Aezg(1)*obj.Aezg(obj.gage.no);%mean(mean(tageswert,2));  % Mittel aller Jahre (2005-2014)
            obj.gage.MQ_musk = MQ(obj.gage.no)*2;%MQ_musk/obj.Aezg(1)*obj.Aezg(obj.gage.no);%mean(mean(tageswert,2));  % Mittel aller Jahre (2005-2014)
            %             obj.gage.MQ = obj.gage.MQ_musk;
            obj.gage.HHQ = HHQ(obj.gage.no)*1.1;%max(max(tageswert));   % Maximaler Abfluss aller Jahre (2005-2014)
            obj.gage.BHQ = HHQ(obj.gage.no)*1.2;
            
        end
        
        function pegel_addAbfluss(obj,gage_no,qgko_abs,qgko_date,qges_abs,gwin_abs,richards,qdir_abs,qifl_abs,qbas_abs)
            %%% Felder leeren
            obj.inflow_wasim = [];
%             obj.discharge_generation = [];
            obj.discharge_qdir = [];%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
            obj.discharge_qifl = [];%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
            obj.discharge_qbas = [];%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];

            obj.gw_infiltration = [];
            obj.inflow{1} = [];
            obj.outflow = [];
            obj.discharge_routed = [];
            obj.routing.fluss = [];
            obj.routing.Qaus = [];
            if isfield(obj.routing, 'Qout') ==1
                obj.routing = rmfield(obj.routing,'Qout');
            end
            
            %%% discharge hydrographs
            % routed discharge (output from WaSiM):
            obj.inflow_wasim(:,1) = qgko_date(:,1);%[repmat(datenum(1990,1,1),24,1);qgko_date(:,1)];
            obj.inflow_wasim(:,2) = qgko_abs(:,gage_no);%[repmat(qgko_abs(1,gage_no),24,1);qgko_abs(:,gage_no)];    % wasim-qgko-time series (before existing lake)
            obj.inflow_wasim(:,3) = obj.inflow_wasim(:,2);          % wasim-qgko-time series (after existing lake)
            % generated discharge (output from WaSiM):
%             obj.discharge_generation(:,1) = qgko_date(:,1);%[repmat(datenum(1990,1,1),24,1);qgko_date(:,1)];
%             obj.discharge_generation(:,2) = qges_abs(:,gage_no);%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
            obj.discharge_qdir(:,2) = qdir_abs(:,gage_no);%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
            obj.discharge_qifl(:,2) = qifl_abs(:,gage_no);%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
            obj.discharge_qbas(:,2) = qbas_abs(:,gage_no);%[repmat(qges_abs(1,gage_no),24,1);qges_abs(:,gage_no)];
                        
            % infiltration from rivers into groundwater:
            obj.gw_infiltration(:,1) = qgko_date(:,1);%[repmat(datenum(1990,1,1),24,1);qgko_date(:,1)];
            obj.gw_infiltration(:,2) = gwin_abs(:,gage_no);%[repmat(qges_abs(1,gage_no),24,1);gwin_abs(:,gage_no)];
            % initial conditions for all gages (no basin + no basin upstream)
            obj.inflow{1} = obj.inflow_wasim(:,1:2);      % no upstream gage -> no routing required
            obj.outflow{1} = obj.inflow{1};           % no basin at the gage -> no basin retention
            obj.discharge_routed(:,1) = obj.inflow{1}(:,1);    % time
%             obj.discharge_routed(:,2) = obj.inflow{1}(:,2)-obj.discharge_generation(:,2);  % routed discharge
            
            obj.discharge_routed(:,2) = obj.inflow{1}(:,2)-obj.discharge_qdir(:,2)-obj.discharge_qifl(:,2)-obj.discharge_qbas(:,2);  % routed discharge
            
            % Richards-Datei
            obj.routing.fluss = richards(1,gage_no).fluss;
            obj.routing.Qaus = richards(1,gage_no).Qaus;
            
        end
        
        function nachfolger(obj)
            %-------------------------------------------------------------------------%
            % function to define the downstream neigbor of a flowpoint (= Nachfolger) %
            %-------------------------------------------------------------------------%
            
            riv_no = [obj.river.number] == obj.gage.river;      % row of the river number
            riv_pegel = sortrows(obj.river(riv_no).gages,5);    % gages of the river
            row = find(riv_pegel(:,1) == obj.gage.no);          % row of the object gage
            
            if obj.gage.no == 1
                % the basin outlet has no downstream neighbor
                obj.neighbors.nachfolger = [];
            else
                if row > 1
                    % if there is a downstream neighbor within the river section
                    obj.neighbors.nachfolger = riv_pegel(row-1,1);
                else
                    % find downstream river section
                    nachfolger=[];
                    while isempty(nachfolger)==1        % (river number changes in every loop)
                        riv_flowpoints = sortrows(obj.river(riv_no).flowpoints,5);	% sort all flow points of the river section with respect to their flow time (ascending order)
                        knot = riv_flowpoints(1,2:3);                               % 1st point: point most downstream
                        flowdir = obj.grids.FLK(knot(1,2),knot(1,1));           % flow direction
                        if flowdir==1; knot2=[knot(1,1)+1,knot(1,2)-1];         % right-up
                        elseif flowdir==2; knot2=[knot(1,1)+1,knot(1,2)];       % right
                        elseif flowdir==4; knot2=[knot(1,1)+1,knot(1,2)+1];     % right-down
                        elseif flowdir==8; knot2=[knot(1,1),knot(1,2)+1];       % down
                        elseif flowdir==16; knot2=[knot(1,1)-1,knot(1,2)+1];    % left-down
                        elseif flowdir==32; knot2=[knot(1,1)-1,knot(1,2)];      % left
                        elseif flowdir==64; knot2=[knot(1,1)-1,knot(1,2)-1];    % left-up
                        elseif flowdir==128; knot2=[knot(1,1),knot(1,2)-1];     % up
                        else disp('error!!')
                        end
                        % next river section:
                        nachfolger_pot = obj.riverpoints(obj.riverpoints(:,2)==knot2(1,1),:);
                        riv_no = nachfolger_pot(nachfolger_pot(:,3)==knot2(1,2),4);
                        % all gages of the new river section:
                        riv_pegel2 = obj.river(riv_no).gages;
                        if isnan(riv_pegel2)==0     % if river section has at least one gage
                            % sort all gages of the river section with respect to their flow time (ascending order)
                            riv_pegel2 = sortrows(obj.river(riv_no).gages,5);
                            % downstream neigbor = most upstream gage of the new river section
                            obj.neighbors.nachfolger = riv_pegel2(end,1);
                            nachfolger = riv_pegel2(end,1);
                        end
                        clear river_depth riv_pegel2
                    end
                end
            end
        end
        
        function pegel = vorgaenger(pegel,i)
            %-------------------------------------------------------------------------%
            % function to define the upstream neigbor(s) of a flowpoint (= Vorgänger) %
            %-------------------------------------------------------------------------%
            
            danach = pegel(1, i).neighbors.nachfolger;
            if isempty(danach)==0
                % if gage i has a downstream neighbor, i is added to the
                % upstream neighbors of that gage
                pegel(1, danach).neighbors.vorgaenger = [pegel(1, danach).neighbors.vorgaenger,pegel(1, i).gage.no];
                pegel(1, danach).pegel_vor_yn = 1;
            end
        end
        
        function pegel = alle_nachfolger(pegel,i)
            %-----------------------------------------------------------------------------%
            % function to define all upstream neigbors of a flowpoint (= alle Nachfolger) %
            %-----------------------------------------------------------------------------%
            
            next_nachfolger = i;
            % loop from one downstream neighbor to the next until gage 1 is
            % reached (all gages on the way are stored in the right order)
            alle_nachfolger = pegel(1, i).neighbors.alle_nachfolger;
            while next_nachfolger~=1
                next_nachfolger = pegel(1, next_nachfolger).neighbors.nachfolger;
                alle_nachfolger = [alle_nachfolger,next_nachfolger];
            end
            pegel(1, i).neighbors.alle_nachfolger = alle_nachfolger;
            
        end
        
        function pegel = becken_ausschluss(pegel,i,gage_outlet,ausschluss)%,method,D_min,D_max)
            %------------------------------------------------%
            % function to define where basins can be located %
            %------------------------------------------------%
            
            % basin can be positioned at a gage, if it is upstream of a
            % specific gage (gage_outlet)
            if any(gage_outlet == pegel(1, i).neighbors.alle_nachfolger)==1
                pegel(1, i).becken_moeglich = 1;
            else
                pegel(1, i).becken_moeglich = 0;
            end
            
            % no basins at positions of lakes (9,10,11,12), extractions (13,14)
            % or existing gages donwstream of the area under investigation
            if any(pegel(1, i).gage.no == ausschluss)==1
                pegel(1, i).becken_moeglich = 0;
            end
            
            
        end
        
        function pegel = becken_vorgaenger(pegel,i)
            %-------------------------------------------------%
            % function to state all upstream basins of a gage %
            %-------------------------------------------------%
            
            % the function is run for every gage with a basin:
            % all gages downstream that gage have the basin upstream
            % -> find all downstream neighbors:
            nachfolger_von_becken = pegel(1, i).neighbors.alle_nachfolger;
            
            for j = nachfolger_von_becken
                % add the basin of the gage i to the upstream basins of that gage
                pegel(1, j).neighbors.becken_vorgaenger = [pegel(1, j).neighbors.becken_vorgaenger, pegel(1, i).gage.no];
                pegel(1, j).becken_vor_yn = 1;  % definition that the gage has a basin upstream (yes = 1)
                pegel(1, j).inflow = [];        % delete inflow (needs to be routed from the upstream basin)
                pegel(1, j).outflow = [];       % delete outflow (= new inflow or discharge after basin retention)
            end
            pegel(1, i).becken_yn = 1;      % definition that the gage has a basin (yes = 1)
            pegel(1, i).outflow = [];       % delete outflow (= new discharge after basin retention)
        end
        
        function abstraction(obj)
            %---------------------------%
            % function for abstractions %
            %---------------------------%
            
            index = (obj.outflow(:,2)>obj.abstr_par.min_abstr).*(obj.outflow(:,2)<=obj.abstr_par.min_abstr+obj.abstr_par.max_abstr/obj.abstr_par.frac_abstr)==1;
            obj.outflow(index,2) = obj.outflow(index,2)-(obj.outflow(index,2)-obj.abstr_par.min_abstr)*obj.abstr_par.frac_abstr;
            index2 = obj.outflow(:,2)>obj.abstr_par.min_abstr+obj.abstr_par.max_abstr/obj.abstr_par.frac_abstr;
            obj.outflow(index2,2) = obj.outflow(index2,2)- obj.abstr_par.max_abstr;
        end
        
        function beckenparametrisierung(obj,k)
            %----------------------------------------------------------------%
            % function to set default values to non-existent basin variables %
            %----------------------------------------------------------------%
            
            %%% basin geometry:
            % basin volume (< 50'000 m^3):
            if isfield(obj.beckenparameter(k),'S_max')==0; obj.beckenparameter(k).S_max = 30000; end
            % storage height (< 4 m):
            if isfield(obj.beckenparameter(k),'h_W')==0; obj.beckenparameter(k).h_W = 3; end
            % form factor (S = a*h^m)
            %             if isfield(obj.beckenparameter(k),'m')==0; obj.beckenparameter.m = 2.5;  end
            %             obj.beckenparameter(k).a = obj.beckenparameter(k).S_max/obj.beckenparameter(k).h_W^obj.beckenparameter(k).m;
            
            
            %%% Berechnung der Öffnungsgrößen (Methode 1)
            obj.beckenparameter(k).AL_min = obj.beckenparameter(k).D_min^2;
            obj.beckenparameter(k).AL_max = obj.beckenparameter(k).D_max^2;
            
            obj.gage.AL_min(k) = obj.beckenparameter(k).AL_min;
            obj.gage.AL_max(k) = obj.beckenparameter(k).AL_max;
            
            %%% Betriebsauslass
            % contraction coefficient (opening):
            %  -> mu_A is defined during the volume-discharge-relation section
            % minimale Öffnungsgröße (to guarantee a baseflow at all time):
            %             if length(obj.gage.AL_min)>1
            %                 obj.beckenparameter(k).AL_min = obj.gage.AL_min(k);
            %                 obj.beckenparameter(k).AL_max = obj.gage.AL_max(k);
            %             else
            %                 obj.beckenparameter(k).AL_min = obj.gage.AL_min;
            %                 obj.beckenparameter(k).AL_max = obj.gage.AL_max;
            %             end
            %             obj.gage.AL_min = [];
            %             obj.gage.AL_max = [];
            %             clear obj.gage.AL_min obj.gage.AL_max
            
            %%% Hochwasserentlastung
            % weir runoff coefficient (weir):
            if isfield(obj.beckenparameter(k),'mu_W')==0; obj.beckenparameter(k).mu_W = 0.72;
            elseif isempty(obj.beckenparameter(k).mu_W)==1; obj.beckenparameter(k).mu_W = 0.72; end
            % freeboard [DIN 19700-12(2004)]:
            if isfield(obj.beckenparameter(k),'f')==0; obj.beckenparameter(k).f = 0.5;
            elseif isempty(obj.beckenparameter(k).f)==1; obj.beckenparameter(k).f = 0.5; end
            % overflow height:
            h_ue = .5;
            if isfield(obj.beckenparameter(k),'h_HWE')==0; obj.beckenparameter(k).h_HWE = obj.beckenparameter(k).h_W + h_ue; end
            
            % weir width:
            max_Q = max(max(obj.gage.BHQ),max(1.2*obj.inflow_wasim(:,3)));
            if isfield(obj.beckenparameter(k),'b_W')==0; obj.beckenparameter(k).b_W = max_Q/(2/3*obj.beckenparameter(k).mu_W*sqrt(2*9.81)*h_ue^(3/2));
            elseif isempty(obj.beckenparameter(k).b_W)==1; obj.beckenparameter(k).b_W = max_Q/(2/3*obj.beckenparameter(k).mu_W*sqrt(2*9.81)*h_ue^(3/2)); end
            
            if obj.beckenparameter(k).b_W < 2
                obj.beckenparameter(k).b_W = 2;
            end
            
            % dam height:
            if isfield(obj.beckenparameter(k),'h_D')==0; obj.beckenparameter(k).h_D = obj.beckenparameter(k).h_HWE + obj.beckenparameter(k).f;
            elseif isempty(obj.beckenparameter(k).h_D)==1; obj.beckenparameter(k).h_D = obj.beckenparameter(k).h_HWE + obj.beckenparameter(k).f;
            end
            
        end
        
        function volumen_abfluss_beziehung(obj,method)
            %---------------------------------------------------------------%
            % function to define the volume-discharge relation of the basin %
            %---------------------------------------------------------------%
            AL = obj.beckenparameter.A_L;
            dist = .1;
            S = dist;
            max_dist = obj.beckenparameter.S_max/100;
            Q = 0;
            
            if method ==1
                % discharge limits:
                % maximal open channel flow = flow under pressure at bankfull stage (h = 0)
                obj.beckenparameter.mu_A = 0.5686;
                obj.beckenparameter.Q_maxF = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*sqrt(AL));
                % maximal flow under pressure = flow under pressure at h = h_W
                if (sqrt(AL)/(obj.beckenparameter.h_W+sqrt(AL)))<0.75
                    psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(obj.beckenparameter.h_W+sqrt(AL)))^2));
                    obj.beckenparameter.mu_A = psi/(sqrt(1+psi*sqrt(AL)/(obj.beckenparameter.h_W+sqrt(AL))));
                else
                    obj.beckenparameter.mu_A = 0.5686;
                end
                obj.beckenparameter.Q_maxL = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*(obj.beckenparameter.h_W+sqrt(AL)));
                % maximal discharge (including emergency spillway)
                obj.beckenparameter.Q_maxG = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*(obj.beckenparameter.h_D-obj.beckenparameter.f+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*obj.beckenparameter.mu_W*obj.beckenparameter.b_W*sqrt(2*9.81)*(obj.beckenparameter.h_D-obj.beckenparameter.f-obj.beckenparameter.h_W)^(3/2);
                
                % open channel flow:
                %                 ShQ = zeros(2,3);
                ShQ(1,:) = [0,0,0];
                ShQ(2,:) = [0,0,obj.beckenparameter.Q_maxF];
                
                % flow under pressure + weir:
                while Q <= obj.gage.HHQ || S <= 2*obj.beckenparameter.S_max;
                    h = (S/obj.beckenparameter.a)^(1/obj.beckenparameter.m);
                    if (sqrt(AL)/(h+sqrt(AL)))<0.75
                        psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h+sqrt(AL)))^2));
                        obj.beckenparameter.mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h+sqrt(AL))));
                    else
                        obj.beckenparameter.mu_A = 0.5686;
                    end
                    if h <= obj.beckenparameter.h_W                             % nur Durckabfluss
                        Q = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)));
                    elseif h <= obj.beckenparameter.h_D-obj.beckenparameter.f   % HWE: Druckabfluss + Wehrüberfall
                        Q = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*obj.beckenparameter.mu_W*obj.beckenparameter.b_W*sqrt(2*9.81)*(h-obj.beckenparameter.h_W)^(3/2);
                    else
                        obj.beckenparameter.ShQ = ShQ;
                        return
                    end
                    ShQ(end+1,:) = [S,h,Q];
                    
                    % select next storage volume for the V-Q-relation
                    dist = min(dist*3/2,max_dist);
                    S = S+dist;
                end
                obj.beckenparameter.ShQ = ShQ;
                
            elseif method ==2
                
                obj.beckenparameter.mu_A = 0.62;
                
                % discharge limits:
                % maximal open channel flow = flow under pressure at bankfull stage (h = 0)
                obj.beckenparameter.Q_maxF = 0;
                % maximal flow under pressure = flow under pressure at h = h_W
                obj.beckenparameter.Q_maxL = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*obj.beckenparameter.h_W);
                % maximal discharge (including emergency spillway)
                obj.beckenparameter.Q_maxG = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*(obj.beckenparameter.h_D-obj.beckenparameter.f))... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*obj.beckenparameter.mu_W*obj.beckenparameter.b_W*sqrt(2*9.81)*(obj.beckenparameter.h_D-obj.beckenparameter.f-obj.beckenparameter.h_W)^(3/2);
                
                % open channel flow:
                ShQ(1,:) = [0,0,0];
                
                % flow under pressure + weir:
                
                while Q <= obj.gage.HHQ || S <= 2*obj.beckenparameter.S_max;
                    h = (S/obj.beckenparameter.a)^(1/obj.beckenparameter.m);
                    if h <= obj.beckenparameter.h_W                             % nur Durckabfluss
                        Q = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*h);
                    elseif h <= obj.beckenparameter.h_D-obj.beckenparameter.f   % HWE: Druckabfluss + Wehrüberfall
                        Q = obj.beckenparameter.mu_A*AL*sqrt(2*9.81*h)... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*obj.beckenparameter.mu_W*obj.beckenparameter.b_W*sqrt(2*9.81)*(h-obj.beckenparameter.h_W)^(3/2);
                    else
                        obj.beckenparameter.ShQ = ShQ;
                        return
                    end
                    ShQ(end+1,:) = [S,h,Q];
                    
                    % select next storage volume for the V-Q-relation
                    dist = min(dist*3/2,max_dist);
                    S = S+dist;
                end
                obj.beckenparameter.ShQ = ShQ;
                
            end
        end
        
        function volumen_abfluss_beziehung_real(obj,method,k)
            %---------------------------------------------------------------%
            % function to define the volume-discharge relation of the basin %
            %---------------------------------------------------------------%
            AL = obj.beckenparameter(k).A_L;
            
            %             S = dist;
            
            Q = 0;
            
            if method ==1
                % discharge limits:
                % maximal open channel flow = flow under pressure at bankfull stage (h = 0)
                mu_A = 0.5686;
                obj.beckenparameter(k).Q_maxF = mu_A*AL*sqrt(2*9.81*sqrt(AL));
                
                % maximal flow under pressure = flow under pressure at h = h_W
                if (sqrt(AL)/(obj.beckenparameter(k).h_W+sqrt(AL)))<0.75
                    psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(obj.beckenparameter(k).h_W+sqrt(AL)))^2));
                    mu_A = psi/(sqrt(1+psi*sqrt(AL)/(obj.beckenparameter(k).h_W+sqrt(AL))));
                else
                    mu_A = 0.5686;
                end
                obj.beckenparameter(k).Q_maxL = mu_A*AL*sqrt(2*9.81*(obj.beckenparameter(k).h_W+sqrt(AL)));
                
                % maximal discharge (including emergency spillway)
                obj.beckenparameter(k).Q_maxG = mu_A*AL*sqrt(2*9.81*(obj.beckenparameter(k).h_D-obj.beckenparameter(k).f+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*obj.beckenparameter(k).mu_W*obj.beckenparameter(k).b_W*sqrt(2*9.81)*(obj.beckenparameter(k).h_D-obj.beckenparameter(k).f-obj.beckenparameter(k).h_W)^(3/2);
                
                % open channel flow:
                ShQ(1,:) = [0,0,0];
                
                ShAs = obj.beckenparameter(k).ShAs;
                
                mu_W = obj.beckenparameter(k).mu_W;
                b_W = obj.beckenparameter(k).b_W;
                h_W = obj.beckenparameter(k).h_W;
                h_D = obj.beckenparameter(k).h_D;
                
                V_W = obj.beckenparameter(k).V_W;
                V_HWE = obj.beckenparameter(k).V_HWE;
                V_D = obj.beckenparameter(k).V_D;
                
                f = obj.beckenparameter(k).f;
                
                last_hW = find(ShAs(:,2)<h_W,1,'last');
                last_hHWE = find(ShAs(:,2)<h_D-f,1,'last');
                last_hD = find(ShAs(:,2)<h_D,1,'last');
                last_hW_check=0;
                last_hHWE_check=0;
                last_hD_check=0;
                
                % flow under pressure + weir:
                i_k=1;
                zeile = 0;
                for z = 1:size(ShAs,1)+sum([isempty(last_hW)==0,isempty(last_hHWE)==0,isempty(last_hD)==0])
                    i_k=i_k+1;
                    zeile = zeile+1;
                    
                    if zeile-last_hW==1 && last_hW_check==0
                        h = h_W;
                        S = V_W;
                        last_hW_check = 1;
                        zeile = zeile-1;
                    elseif zeile>last_hHWE && last_hHWE_check==0
                        h = h_D-f;
                        S = V_HWE;
                        last_hHWE_check = 1;
                        zeile = zeile-1;
                    elseif zeile>last_hD && last_hD_check==0
                        h = h_D;
                        S = ShAs(find(abs(ShAs(:,2)-h)<1e-10,1,'first'),4);
                        last_hD_check = 1;
                        zeile = zeile-1;
                    else
                        
                        h = ShAs(zeile,2); % erste Zeile: h=0
                        S = ShAs(zeile,4);
                    end
                    if (sqrt(AL)/(h+sqrt(AL)))<0.75
                        psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h+sqrt(AL)))^2));
                        mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h+sqrt(AL))));
                    else
                        mu_A = 0.5686;
                    end
                    if h <= h_W % nur Druckabfluss
                        Q = mu_A*AL*sqrt(2*9.81*(h+sqrt(AL))); % bei h = 0: maximaler Freispiegelabfluss (für zeile=1)
                        %         bereich = 1;
                    elseif h <= h_D-f   % HWE: Druckabfluss + Wehrüberfall
                        Q = mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                        %         bereich = 2;
                    elseif h<=h_D
                        Q = mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                        %         bereich = 3;
                    else
                        Q = mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                        %         bereich = 4;
                        ShQ(i_k,1) = S-ShAs(1,4);
                        ShQ(i_k,2) = h;
                        ShQ(i_k,3) = Q;
                        obj.beckenparameter(k).ShQ = ShQ;
                        continue
                    end
                    
                    ShQ(i_k,1) = S-ShAs(1,4);
                    ShQ(i_k,2) = h;
                    ShQ(i_k,3) = Q;
                    %                     [ShAs(zeile,4),h,Q];
                    
                end
                
                h=h+3;
                S = ShQ(end,1)+3*ShAs(end,3);
                mu_A = 0.5686;
                Q = mu_A*AL*sqrt(2*9.81*(h+sqrt(AL)))... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                ShQ(end+1,:) = [S,h,Q];
                obj.beckenparameter(k).ShQ = ShQ;
                
                
                
                
            elseif method ==2
                
                mu_A = 0.62;
                
                mu_W = obj.beckenparameter(k).mu_W;
                b_W = obj.beckenparameter(k).b_W;
                h_W = obj.beckenparameter(k).h_W;
                h_D = obj.beckenparameter(k).h_D;
                f = obj.beckenparameter(k).f;
                
                % discharge limits:
                % maximal open channel flow = flow under pressure at bankfull stage (h = 0)
                obj.beckenparameter(k).Q_maxF = 0;
                % maximal flow under pressure = flow under pressure at h = h_W
                obj.beckenparameter(k).Q_maxL = mu_A*AL*sqrt(2*9.81*h_W);
                % maximal discharge (including emergency spillway)
                obj.beckenparameter(k).Q_maxG = mu_A*AL*sqrt(2*9.81*(h_D-f))... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*mu_W*b_W*sqrt(2*9.81)*(h_D-f-h_W)^(3/2);
                
                % open channel flow:
                %                 ShQ(1,:) = [0,0,0];
                
                ShAs = obj.beckenparameter(k).ShAs;
                
                
                
                % flow under pressure + weir:
                
                for zeile = 1:size(obj.beckenparameter(k).ShAs,1)
                    %                 while Q <= obj.gage.HHQ || S <= 2*obj.beckenparameter(k).S_max;
                    h = ShAs(zeile,2);
                    %                     h = obj.beckenparameter(k).ShAs(zeile,2);%(S/obj.beckenparameter(k).a)^(1/obj.beckenparameter(k).m);
                    if h <= h_W                             % nur Durckabfluss
                        Q = mu_A*AL*sqrt(2*9.81*h);
                    elseif h <= h_D-f   % HWE: Druckabfluss + Wehrüberfall
                        Q = mu_A*AL*sqrt(2*9.81*h)... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                    else
                        obj.beckenparameter(k).ShQ = ShQ;
                        return
                    end
                    ShQ(zeile,:) = [ShAs(zeile,4)-ShAs(1,4),h,Q];
                    
                    % select next storage volume for the V-Q-relation
                    %                     dist = min(dist*3/2,max_dist);
                    %                     S = S+dist;
                end
                obj.beckenparameter(k).mu_A = mu_A;
                
                
                dist = (ShQ(end,1)-ShQ(end-1,1))/2;
                max_dist = obj.beckenparameter(k).S_max/50;
                
                
                % Hochwasserentlastung
                while ShQ(end,3) <= max(obj.gage.HHQ,1.5*max(obj.inflow_wasim(:,3))) && ShQ(end,1) <= 2.5*obj.beckenparameter(k).S_max
                    % select next storage volume for the V-Q-relation
                    dist = min(dist*3/2,max_dist);
                    S = ShQ(end,1)+dist;
                    h = h+dist/ShAs(end,3);
                    
                    
                    if h <= h_W                             % nur Durckabfluss
                        Q = mu_A*AL*sqrt(2*9.81*h);
                    elseif h <= h_D-f   % HWE: Druckabfluss + Wehrüberfall
                        Q = mu_A*AL*sqrt(2*9.81*h)... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                    else
                        h=h+3;
                        S = ShQ(end,1)+3*ShAs(end,3);
                        Q = mu_A*AL*sqrt(2*9.81*h)... % Druckabfluss bei maximalem Wasserstand
                            + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                        ShQ(end+1,:) = [S,h,Q];
                        obj.beckenparameter(k).ShQ = ShQ;
                        return
                    end
                    
                    %                         S = ShQ(end,1) + ShAs(end,3)*.05;
                    ShQ(end+1,:) = [S,h,Q];
                    
                end
                
                % Extremabfluss, damit die Tabelle in jedem Fall
                % ausreicht (wird nie sinnvoll eintreten --> Becken
                % fällt dann eh raus)
                h=h+3;
                S = ShQ(end,1)+3*ShAs(end,3);
                Q = mu_A*AL*sqrt(2*9.81*h)... % Druckabfluss bei maximalem Wasserstand
                    + 2/3*mu_W*b_W*sqrt(2*9.81)*(h-h_W)^(3/2);
                ShQ(end+1,:) = [S,h,Q];
                
                
                obj.beckenparameter(k).ShQ = ShQ;
                
                
                
            end
        end
        
        function beckenoutput_wasim(obj)
            %--------------------------------------------------------%
            % function for data table output (for import into WaSiM) %
            %--------------------------------------------------------%
            
            %             step = (obj.gage.BHQ-obj.beckenparameter.Q_maxF)/25;
            anfang = find(obj.beckenparameter.ShQ(:,3)>=obj.beckenparameter.Q_maxF,1,'first');
            becken_wasim(1,1) = obj.beckenparameter.ShQ(anfang,1);
            becken_wasim(1,2) = obj.beckenparameter.ShQ(anfang,3);
            %             for ii = anfang+1:length(obj.beckenparameter.ShQ(:,3))
            %                 if obj.beckenparameter.ShQ(ii,3) < 2*obj.gage.BHQ
            %                     if (obj.beckenparameter.ShQ(ii,3)-becken_wasim(end,2))>step
            %                         becken_wasim(end+1,1) = obj.beckenparameter.ShQ(ii,1);
            %                         becken_wasim(end,2) = obj.beckenparameter.ShQ(ii,3);
            %                     end
            %                 else
            %                     break
            %                 end
            %             end
            step = 1;
            while obj.beckenparameter.ShQ(step+anfang-1,3)<2*obj.gage.BHQ
                becken_wasim(end+1,1) = obj.beckenparameter.ShQ(step+anfang,1);
                becken_wasim(end,2) = obj.beckenparameter.ShQ(step+anfang,3);
                step = min(ceil(7/6*step),step+30);
                if step>length(obj.beckenparameter.ShQ)-anfang
                    break
                end
            end
            becken_wasim(end+1,1) = obj.beckenparameter.ShQ(end,1);
            becken_wasim(end,2) = 2*obj.gage.HHQ;
            obj.beckenparameter.wasim = becken_wasim;
        end
        
    end
end
