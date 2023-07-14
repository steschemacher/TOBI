classdef Flusspunkt<handle
    
    properties
        ID
        coordinates
        nachfolger
        vorgaenger
        grids
        flussnetz
        elevation
        fliesszeit
    end
    
    methods
        function obj = Flusspunkt(i,grids,riv_point)
            obj.grids.FLK = grids.FLK;          % flow directions
            obj.grids.RIV = grids.RIV;          % Flussnetz (Matrix)
            obj.flussnetz = riv_point;          % ID + Koordinaten
            obj.coordinates = obj.flussnetz(i,2:3); % Objekt-Koordinaten
            obj.elevation = grids.DGM(obj.coordinates(1,2),obj.coordinates(1,1));
            obj.fliesszeit = grids.FZS(obj.coordinates(1,2),obj.coordinates(1,1));
            obj.ID = i;
            obj.vorgaenger = [];
        end
        function fun_nachfolger(obj)
            % function to define the downstream neigbor of a flowpoint
            knot = obj.coordinates;
            flowdir = obj.grids.FLK(knot(1,2),knot(1,1));           % Richtung, in die es an diesem Punkt fließt
            if flowdir==1; knot2=[knot(1,1)+1,knot(1,2)-1];         % rechts-oben
            elseif flowdir==2; knot2=[knot(1,1)+1,knot(1,2)];       % rechts
            elseif flowdir==4; knot2=[knot(1,1)+1,knot(1,2)+1];     % rechts-unten
            elseif flowdir==8; knot2=[knot(1,1),knot(1,2)+1];       % unten
            elseif flowdir==16; knot2=[knot(1,1)-1,knot(1,2)+1];    % links-unten
            elseif flowdir==32; knot2=[knot(1,1)-1,knot(1,2)];      % links
            elseif flowdir==64; knot2=[knot(1,1)-1,knot(1,2)-1];    % links-oben
            elseif flowdir==128; knot2=[knot(1,1),knot(1,2)-1];     % oben
            else
                disp('error!!')
            end
            nachfolger_pot = obj.flussnetz(obj.flussnetz(:,2)==knot2(1,1),:);
            obj.nachfolger = nachfolger_pot(nachfolger_pot(:,3)==knot2(1,2),1);
        end  
        function fun_vorgaenger(flusspunkte,i)
            % function to define the upstream neigbor(s) of a flowpoint
            danach = flusspunkte(1, i).nachfolger;
            if isempty(danach)==0;
                flusspunkte(1, danach).vorgaenger = [flusspunkte(1, danach).vorgaenger,flusspunkte(1, i).ID];
            end
        end
    end
end