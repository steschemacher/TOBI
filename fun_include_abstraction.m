function pegel = fun_include_abstraction(b_nr,abstraction,min_abstr,max_abstr,frac_abstr,pegel)

% function to include abstractions to the routing structure

pegel(b_nr).abstraction_yn = 1;

if pegel(b_nr).becken_vor_yn==1;
pegel(b_nr).inflow_wasim(:,2) = pegel(b_nr).inflow_wasim(:,3)+abstraction.abstraction;
end
pegel(b_nr).abstr_par.min_abstr = min_abstr;
pegel(b_nr).abstr_par.max_abstr = max_abstr;
pegel(b_nr).abstr_par.frac_abstr = frac_abstr;
