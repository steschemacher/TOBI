function [zeitraum_B,zeitraumA_B] = fun_zeitraum_B(abfluss,QmaxF)

zeitraum_B = find(abfluss>QmaxF,1):find(abfluss(find(abfluss==max(abfluss),1,'first'):end)>QmaxF,1,'last')+find(abfluss==max(abfluss),1,'first')-1;
zeitraumA_B = find(abfluss>QmaxF,1):find(abfluss==max(abfluss),1,'first');