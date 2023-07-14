function [pegel,max_gage] = f2_1__define_gage_properties(grids,river,A_ezg,river_points,Q_mean,HHQ_Ereignisse,routing)

max_gage = max(max(grids.PUR));

for i = 1:max_gage
    % Allgemeine Pegeleigenschaften zuweisen
    pegel(i) = Pegel(i,grids,river,A_ezg,river_points,Q_mean,HHQ_Ereignisse);
    
    % Nachfolger bestimmen
    pegel(i).nachfolger;
end

for i = 1:max_gage
    % find the upstream neighbor(s) of every gage
    pegel = vorgaenger(pegel,i);
    % find all downstream neigbors of every gage (to the basin outlet)
    pegel = alle_nachfolger(pegel,i); % includes setting possible gages only upstream of the gage "regional"
    
    if i>1
        pegel(i).routing = routing(i);
    end
    pegel(i).neighbors.alle_vorgaenger = [];
end
for i = 1:max_gage
    for j = pegel(i).neighbors.alle_nachfolger
        pegel(j).neighbors.alle_vorgaenger = [pegel(j).neighbors.alle_vorgaenger,i];
    end
end

