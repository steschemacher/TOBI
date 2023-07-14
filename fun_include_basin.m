function pegel = fun_include_basin(b_nr,basin,pegel)

% function to include basins to the routings structure

pegel(b_nr).becken_yn = 1;
pegel(b_nr).VQ_fest_yn = 1;
pegel(b_nr).beckenparameter.Q_maxF = 0;

pegel(b_nr).beckenparameter.ShQ(:,1) = basin.ShQ(:,1);
pegel(b_nr).beckenparameter.ShQ(:,3) = basin.ShQ(:,3);
 
pegel(b_nr).inflow_wasim(:,2) = basin.basin_rule(:,10);
% if b_nr ==9
%     pegel(b_nr).beckenparameter.S0 = 320E06;
% else
    pegel(b_nr).beckenparameter.S0 = basin.S0;
% end

% external influences on the lake volume: -Irrig + QGW - ETR + PREC - Abstr + Err
% if b_nr<11
pegel(b_nr).beckenparameter.difference = (-basin.basin_rule(:,16)+basin.basin_rule(:,19)...
    -basin.basin_rule(:,22)+basin.basin_rule(:,25)...
    -basin.basin_rule(:,31)+basin.basin_rule(:,34));
% else
%     pegel(b_nr).beckenparameter.difference = basin.basin_rule(:,16);
% end

