%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%               Funktion zum Starten den Beckenszenarien                 %
%                                                                        %
%     Das Programm enthält alle einzugsgebietsspezifischen Parameter     %
%       - Ordnernamen                                                    %
%       - Beckenanzahlen                                                 %
%       - Routingparameter                                               %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f1__catchment_input(Gebiet,TG,...
    E_no_opt,Beckenkombination_speichern,...
    mOpt,opt_parameter,method)

path_allg = 'PFAD CODE';
path_becken = 'PFAD BECKEN LOCASIN';
% method = 1;

%%% Abhängig vom Untersuchungsgebiet
if Gebiet == 1
%     A_ges = 1167.9;   % Einzugsgebietsgröße (immer Gesamtgebiet)
%     pegel_ausschluss = [1,3,9];   % Pegelnummern, an denen keine Becken angeordnet werden dürfen
    beobachtungspunkte = [1:19];    % Pegelpunkte für "regional"-Auswertung
    
elseif Gebiet == 2
%     A_ges = 392.17;   % Einzugsgebietsgröße (immer Gesamtgebiet)
%     pegel_ausschluss = [1,2,6,9,10,11,12];   % Pegelnummern, an denen keine Becken angeordnet werden dürfen
    beobachtungspunkte = [1,2,6,13,296,247,217,232,268,142,79];    % Pegelpunkte für "regional"-Auswertung
    
elseif Gebiet == 3
%     A_ges = 104;   % Einzugsgebietsgröße (immer Gesamtgebiet)
%     pegel_ausschluss = [1,2];   % Pegelnummern, an denen keine Becken angeordnet werden dürfen
    beobachtungspunkte = [1,2,13,296,247,217,232,268,142,79];    % Pegelpunkte für "regional"-Auswertung
    
elseif Gebiet == 4
%     A_ges = 91.1;   % Einzugsgebietsgröße (immer Gesamtgebiet)
%     pegel_ausschluss = 1;   % Pegelnummern, an denen keine Becken angeordnet werden dürfen
    beobachtungspunkte = [1,2,3,19,115,192,295,325];    % Pegelpunkte für "regional"-Auswertung

end



%% Funktion starten
f2__wirksamkeitsanalyse(Gebiet,TG,...
    E_no_opt,Beckenkombination_speichern,...
    mOpt,opt_parameter,...
    path_allg,path_becken,...
    method,beobachtungspunkte)