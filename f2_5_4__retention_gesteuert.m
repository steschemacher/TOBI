function pegel = f2_5_4__retention_gesteuert(pegel,pegel_no,k)
clear A_L AL AL_alt ans becken_outflow h_speicher hcum i inflow mu_A outflow psi A_speicher Qcum Qmax Qsort Qsort2 T_begin T_ende V_speicher Vcum Q_speicher
inflow = pegel(pegel_no).inflow;

Qsort = sort(inflow{k}(:,2),'descend');
Qsort2 = (Qsort)*3600;
Qcum = cumsum(Qsort2);
hcum = 1:length(Qcum);
Vcum = Qcum-Qsort2.*hcum';
V = pegel(pegel_no).beckenparameter.S_max;
line1 = find(Vcum<=V,1,'last');
if line1==length(Qsort)
    Qmax = pegel(pegel_no).gage.MQ*2;
    Qmax = max(Qmax,min(Qsort));
%     keyboard
else
    line2 = find(Vcum>=V,1,'first');
    V1 = Vcum(line1);
    V2 = Vcum(line2);
    Qmax1 = Qsort(line1);
    Qmax2 = Qsort(line2);
    Qmax = (Qmax2-Qmax1)*(V-V1)/(V2-V1)+Qmax1;
    Qmax = max(Qmax,min(Qsort));
end

Qmax = max(Qmax,pegel(pegel_no).gage.MQ*2);

T_begin = find(inflow{k}(:,2)>=Qmax,1,'first');
T_ende = find(inflow{k}(:,2)>=Qmax,1,'last');


% Speichervolumen und Wasserstand
Q_speicher = (inflow{k}(:,2)-Qmax)*3600;
Q_speicher(Q_speicher<0) = 0;
V_speicher = cumsum(Q_speicher);
V_speicher = V_speicher(1:T_ende,1);
h_speicher =  interp1(pegel(pegel_no).beckenparameter.ShAs(:,4),pegel(pegel_no).beckenparameter.ShAs(:,2),V_speicher);

outflow(1:T_begin-1,1) = inflow{k}(1:T_begin-1,2);
outflow(T_begin:T_ende,1) = Qmax;
pegel(pegel_no).beckenparameter(k).Q_maxF = inflow{k}(find(V_speicher>0,1,'first'),2);
pegel(pegel_no).beckenparameter(k).Q_maxL = Qmax;

for i = 1:length(outflow)
    if i>1
        AL_alt = A_L(i-1);
    else
        AL_alt = 2;
    end
    if h_speicher(i) == 0
        AL = 10;
    else
        if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
            psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
            mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
        else
            mu_A = 0.5686;
        end
        AL = outflow(i,1)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
        while abs(AL-AL_alt)>.05
            AL_alt = AL;
            if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
                psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
                mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
            else
                mu_A = 0.5686;
            end
            AL = outflow(i,:)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
            
        end
    end
    A_L(i,1) = AL;
end
i_max = i;
while V_speicher(i)>0 && i<length(inflow{k}(:,2))%i+1:i+floor(V_speicher(i)/((Qmax-25)*3600))
    i=i+1;
    outflow(i,1) = Qmax;
    V_speicher(i,1) = V_speicher(i-1,1)+ inflow{k}(i,2)*3600 - outflow(i,1)*3600;
    h_speicher(i,1) = interp1(pegel(pegel_no).beckenparameter.ShAs(:,4),pegel(pegel_no).beckenparameter.ShAs(:,2),V_speicher(i,1));
    if i>1
        AL_alt = A_L(i-1);
    else
        AL_alt = 2;
    end
    if h_speicher(i) == 0
        AL = 10;
    else
        if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
            psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
            mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
        else
            mu_A = 0.5686;
        end
        AL = outflow(i,1)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
        while abs(AL-AL_alt)>.05
            AL_alt = AL;
            if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
                psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
                mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
            else
                mu_A = 0.5686;
            end
            AL = outflow(i,:)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
            
        end
    end
    A_L(i,1) = AL;
    
    if V_speicher(i,1)<0
        outflow(i,1) = inflow{k}(i,2);
        V_speicher(i,1) = 0;
        h_speicher(i,1) = 0;
        A_L(i,1) = 10;
    end
end
i_diff = i-i_max;
while V_speicher(i)>0 && i<length(inflow{k}(:,2))
  
i = i+1;

outflow(i,1) = inflow{k}(i-i_diff,2);
    V_speicher(i,1) = V_speicher(i-1,1)+ inflow{k}(i,2)*3600 - outflow(i,1)*3600;
    h_speicher(i,1) = interp1(pegel(pegel_no).beckenparameter.ShAs(:,4),pegel(pegel_no).beckenparameter.ShAs(:,2),V_speicher(i,1));
    if i>1
        AL_alt = A_L(i-1);
    else
        AL_alt = 2;
    end
    if h_speicher(i) == 0
        AL = 10;
    else
        if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
            psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
            mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
        else
            mu_A = 0.5686;
        end
        AL = outflow(i,1)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
        while abs(AL-AL_alt)>.05
            AL_alt = AL;
            if (sqrt(AL)/(h_speicher(i)+sqrt(AL)))<0.75
                psi = 1/(1+.64*sqrt(1-(sqrt(AL)/(h_speicher(i)+sqrt(AL)))^2));
                mu_A = psi/(sqrt(1+psi*sqrt(AL)/(h_speicher(i)+sqrt(AL))));
            else
                mu_A = 0.5686;
            end
            AL = outflow(i,:)/(mu_A*sqrt(2*9.81 * (h_speicher(i) + sqrt(AL_alt))));
            
        end
    end
    A_L(i,1) = AL;
    

% if A_L(i,1)<2.9
%     A_L(i,1) = 2.9;
%     % A_L(i,1) = A_L(i-1);
%         becken_outflow = mu_A*A_L(i,1)*sqrt(2*9.81*(h_speicher(i-1)+sqrt(A_L(i,1))));
%         outflow(i,1) = becken_outflow;
%         V_speicher(i,1) = V_speicher(i-1,1)+ inflow{k}(i,2)*3600 - becken_outflow*3600;
%         h_speicher(i,1) = interp1(pegel(pegel_no).beckenparameter.ShAs(:,4),pegel(pegel_no).beckenparameter.ShAs(:,2),V_speicher(i,1));
% end
        
end
if i<length(inflow{k}(:,2))
    outflow(i+1:length(inflow{k}(:,2)),1) = inflow{k}(i+1:length(inflow{k}(:,2)),2);
    V_speicher(i+1:length(inflow{k}(:,2)),1) = 0;
    h_speicher(i+1:length(inflow{k}(:,2)),1) = 0;
     A_L(i+1:length(inflow{k}(:,2)),1) = A_L(i);
end

pegel(pegel_no).outflow{k}(:,1) = inflow{k}(:,1);
pegel(pegel_no).outflow{k}(:,2) = outflow(:,1);
pegel(pegel_no).becken(k).S(:,2) = V_speicher(:,1);
pegel(pegel_no).becken(k).S_max = max(V_speicher(:,1));
pegel(pegel_no).becken(k).h(:,2) = h_speicher(:,1);
pegel(pegel_no).becken(k).A_L(:,2) = A_L(:,1);
pegel(pegel_no).beckenparameter(k).A_L = Qmax;
pegel(pegel_no).beckenfuellung(k) = max(V_speicher)/pegel(pegel_no).beckenparameter(k).S_max;

%%
% figure
% 
% subplot(2,1,1)
% hold on
% plot(inflow{k}(:,2))
% plot(outflow)
% grid on
% grid minor
% set(gca,'xlim',[0,500])
% 
% subplot(2,1,2)
% hold on
% yyaxis left
% % plot(V_speicher)
% plot(h_speicher)
% 
% yyaxis right
% plot(A_L,'.-')
% % plot(h_speicher)
% grid on
% grid minor
% set(gca,'xlim',[0,500])