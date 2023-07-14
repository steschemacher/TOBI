function pegel = fun_S_gesamt(pegel,pegel_no)

% function to sum up the retention volume upstream of a particular gage

S_gesamt_vorgaenger(1:3) = [0,0,0];
if pegel(pegel_no).becken_vor_yn ==1
    for i = pegel(pegel_no).neighbors.vorgaenger
        S_gesamt_vorgaenger(1) = S_gesamt_vorgaenger(1)+pegel(i).gage.S_gesamt(end,1);
        S_gesamt_vorgaenger(2) = S_gesamt_vorgaenger(2)+pegel(i).gage.S_gesamt(end,2);
        S_gesamt_vorgaenger(3) = S_gesamt_vorgaenger(3)+pegel(i).gage.S_gesamt(end,3);
    end
end

if isprop(pegel(pegel_no),'becken')==1
    if isfield(pegel(pegel_no).becken,'S_max')==1 && isfield(pegel(pegel_no).beckenparameter,'S_max')==1
        for k = 1:length(pegel(pegel_no).beckenparameter)
        pegel(pegel_no).gage.S_gesamt(k,1) = S_gesamt_vorgaenger(1)+sum([pegel(pegel_no).becken(k).S_max]); % gesamtes gespeichertes Volumen MIT HWE
        pegel(pegel_no).gage.S_gesamt(k,2) = S_gesamt_vorgaenger(2)+sum(min(vertcat([pegel(pegel_no).becken(k).S_max],[pegel(pegel_no).beckenparameter(k).S_max]))); % gesamtes gespeichertes Volumen OHNE HWE
        pegel(pegel_no).gage.S_gesamt(k,3) = S_gesamt_vorgaenger(3)+sum([pegel(pegel_no).beckenparameter(k).S_max]);
        
        S_gesamt_vorgaenger(1) = pegel(pegel_no).gage.S_gesamt(k,1);
        S_gesamt_vorgaenger(2) = pegel(pegel_no).gage.S_gesamt(k,2);
        S_gesamt_vorgaenger(3) = pegel(pegel_no).gage.S_gesamt(k,3);

        end
     else
        pegel(pegel_no).gage.S_gesamt(end,1:3) = S_gesamt_vorgaenger(1:3);
    end
else
    pegel(pegel_no).gage.S_gesamt(end,1:3) = S_gesamt_vorgaenger(1:3);
end
