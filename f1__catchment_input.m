%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%               Funktion zum Starten den Beckenszenarien                 %
%                                                                        %
%     Das Programm enth�lt alle einzugsgebietsspezifischen Parameter     %
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

%%% Abh�ngig vom Untersuchungsgebiet
if Gebiet == 1
%     A_ges = 1167.9;   % Einzugsgebietsgr��e (immer Gesamtgebiet)
%     pegel_ausschluss = [1,3,9];   % Pegelnummern, an denen keine Becken angeordnet werden d�rfen
    beobachtungspunkte = [1:19];    % Pegelpunkte f�r "regional"-Auswertung
    
elseif Gebiet == 2
%     A_ges = 392.17;   % Einzugsgebietsgr��e (immer Gesamtgebiet)
%     pegel_ausschluss = [1,2,6,9,10,11,12];   % Pegelnummern, an denen keine Becken angeordnet werden d�rfen
    beobachtungspunkte = [1,2,6,13,296,247,217,232,268,142,79];    % Pegelpunkte f�r "regional"-Auswertung
    
elseif Gebiet == 3
%     A_ges = 104;   % Einzugsgebietsgr��e (immer Gesamtgebiet)
%     pegel_ausschluss = [1,2];   % Pegelnummern, an denen keine Becken angeordnet werden d�rfen
    beobachtungspunkte = [1,2,13,296,247,217,232,268,142,79];    % Pegelpunkte f�r "regional"-Auswertung
    
elseif Gebiet == 4
%     A_ges = 91.1;   % Einzugsgebietsgr��e (immer Gesamtgebiet)
%     pegel_ausschluss = 1;   % Pegelnummern, an denen keine Becken angeordnet werden d�rfen
    beobachtungspunkte = [1,2,3,19,115,192,295,325];    % Pegelpunkte f�r "regional"-Auswertung

end



%% Funktion starten
f2__wirksamkeitsanalyse(Gebiet,TG,...
    E_no_opt,Beckenkombination_speichern,...
    mOpt,opt_parameter,...
    path_allg,path_becken,...
    method,beobachtungspunkte)