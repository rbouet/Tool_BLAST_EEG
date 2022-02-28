function Build_Rapport_V02()
%% ____________ Resultats BLAST-EEG ____________

fprintf('\n\n\n\n\n\n\n\n')

%%
% Le protocole BLAST-EEG consiste à réaliser plusieurs sessions du test Bron-Lyon-Attention-Stability-Test couplé et synchronisé à un électroencéphalogramme.
% Il a pour objectif de déterminer si des  caractéristiques électrophysiologiques normales (ex analyse spectrale, ratio théta/béta pour les enfants TDAH) ou pathologiques (ex anomalies EEG intercritiques pour les enfants épileptiques) peuvent être responsable de difficultés attentionnelles.

fprintf('\n\n\n\n\n\n')

global GUI

[~, nameTRC, extTRC] = fileparts(GUI.BLAST_Object.Path_TRC);
fprintf('Nom du fichier EEG : %s\n', [nameTRC, extTRC]);


fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')


%% 
% *Global*


%%
datas = [sum(strcmpi(GUI.BLAST_Object.AIC.label, 'a'));
         sum(strcmpi(GUI.BLAST_Object.AIC.label, 'c'));
         sum(strcmpi(GUI.BLAST_Object.AIC.label, 'g'));
         sum(strcmpi(GUI.BLAST_Object.AIC.label, 'i'));
         sum(strcmpi(GUI.BLAST_Object.AIC.label, 'e'));
         length(GUI.BLAST_Object.AIC.code_label)/2];
T = array2table(datas,...
    'VariableNames',{'Epileptic_Events'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other',...
                 'Total'});
disp(T)
clear datas
fprintf('\n\n');

datas = [GUI.BLAST_Object.session_duration/60;
         GUI.BLAST_Object.session_duration_without_eye/60];
T = array2table(datas,...
    'VariableNames',{'Duration'},...
    'RowNames', {'Total','Without Eyes blocs'});
disp(T)
clear datas



fprintf('\n\nDécoupage des blocs sur les codes: %s', GUI.BLAST_Object.trig_start_bloc)

% Figure Global
time_bloc = [];
ampl_bloc = [];
for xi = 1 : numel(GUI.BLAST_Object.bloc)
    time_bloc = [time_bloc GUI.BLAST_Object.bloc{xi}.start_stop(1) GUI.BLAST_Object.bloc{xi}.start_stop(1) GUI.BLAST_Object.bloc{xi}.start_stop(2) GUI.BLAST_Object.bloc{xi}.start_stop(2)];
    ampl_bloc = [ampl_bloc 0 1 1 0];
end, clear xi


GUI.fig_AIC_Pos.main = figure('Color', GUI.Colors(1,:), 'MenuBar', 'no', 'position',[20 100 1000 500],...
                  'Tag', 'Pos_AIC_fig');
GUI.fig_AIC_Pos.Bloc =  fill(time_bloc, ampl_bloc, 'k-');
GUI.fig_AIC_Pos.Bloc.Parent.Color = GUI.Colors(4,:);
GUI.fig_AIC_Pos.Bloc.Parent.XLabel.String = 'Time (s)';
GUI.fig_AIC_Pos.Bloc.Parent.XLabel.FontSize = 20;
GUI.fig_AIC_Pos.Bloc.Parent.YTick = [];
GUI.fig_AIC_Pos.Bloc.Parent.YLim = [-0.9 1.2];

hold on
% Plot open eyes/ close eyes
if ~isempty(GUI.BLAST_Object.OY.seconde)
    GUI.fig_AIC_Pos.OY = fill([GUI.BLAST_Object.OY.seconde GUI.BLAST_Object.OY.seconde(2:-1:1)], ampl_bloc(2:5), '-', 'FaceColor', [0.5 0.5 0.5]);
end
if ~isempty(GUI.BLAST_Object.CY.seconde)
    GUI.fig_AIC_Pos.CY = fill([GUI.BLAST_Object.CY.seconde GUI.BLAST_Object.CY.seconde(2:-1:1)], ampl_bloc(2:5), '-', 'FaceColor', [0.5 0.5 0.5]);
end

% Plot AIC
time_aic  = GUI.BLAST_Object.AIC.sample(1:2:end);
label_aic = cell2mat(GUI.BLAST_Object.AIC.label);
label_aic = label_aic(1:2:end);
ampl_aic = 0.5*rand(1,length(time_aic));
time_aic_gen = time_aic(find(label_aic == 'a'));
ampl_aic_gen = ampl_aic(find(label_aic == 'a'));
time_aic_foc = time_aic(find(label_aic == 'c'));
ampl_aic_foc = ampl_aic(find(label_aic == 'c'));
time_other   = time_aic(find(label_aic == 'e'));
ampl_other   = ampl_aic(find(label_aic == 'e'));
time_sei_gen = time_aic(find(label_aic == 'g'));
ampl_sei_gen = ampl_aic(find(label_aic == 'g'));
time_sei_foc = time_aic(find(label_aic == 'i'));
ampl_sei_foc = ampl_aic(find(label_aic == 'i'));

GUI.fig_AIC_Pos.AIC = [plot(time_aic_gen(1:2:end), ampl_aic_gen(1:2:end), 'ro'),...
                       plot(time_aic_foc(1:2:end), ampl_aic_foc(1:2:end), 'r*'),...
                       plot(time_other(1:2:end), ampl_other(1:2:end), 'g*'),...
                       plot(time_sei_gen(1:2:end), ampl_sei_gen(1:2:end), 'bo'),...
                       plot(time_sei_foc(1:2:end), ampl_sei_foc(1:2:end), 'b*')];
                   
t_ref = GUI.fig_AIC_Pos.Bloc.XData(1);
text(t_ref-(t_ref*0.6), -0.3, '* AIC General', 'FontSize',14, 'Color', 'r')
text(t_ref-(t_ref*0.6), -0.4, 'o AIC Focal', 'FontSize',14, 'Color', 'r')
text(t_ref-(t_ref*0.6), -0.5, '* Other', 'FontSize',14, 'Color', [0 1 0])
text(t_ref-(t_ref*0.6), -0.6, 'o Seizure General', 'FontSize',14, 'Color', 'b')
text(t_ref-(t_ref*0.6), -0.7, '* Seizure Focal', 'FontSize',14, 'Color', 'b')
clear t_ref

clear time_aic time_aic_foc time_aic_gen time_other time_sei_foc time_sei_gen 
clear ampl_aic ampl_aic_foc ampl_aic_gen ampl_other ampl_sei_foc ampl_sei_gen


% infos text
for xi_bloc = 1 : numel(GUI.BLAST_Object.bloc)
    text(time_bloc(1 + ((xi_bloc-1))*4), -0.3, [num2str(round(time_bloc(3 + ((xi_bloc-1))*4) - time_bloc(1 + ((xi_bloc-1))*4))) ' sec'], 'FontSize',14)
    text(time_bloc(1 + ((xi_bloc-1))*4), -0.5, [num2str(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB)) ' events'], 'FontSize',14)
    text(time_bloc(1 + ((xi_bloc-1))*4), -0.7, [GUI.BLAST_Object.bloc{xi_bloc}.version], 'FontSize',14)
end, clear xi_bloc

snapnow
delete(GUI.fig_AIC_Pos.main)

%%
% *Visualisation de l'ensemble du protocole* :   
%
% *Colonne verte*  : 2 min de repos yeux ouverts;
% *Colonne orange* : 2 min de repos yeux fermés; 
% *Colonne noire* :  test BLAST; 
% *Points* ou *Étoiles rouges*, *verts* ou *bleus* : répartitions des anomalies EEG critiques ou intercritiques
% Durée et nombre d'items des sessions 

fprintf('\n\n\n\n\n\n\n\n');
 





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOC 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% *BLOC 1*

xi_bloc = 1;

str = ['Bloc ', num2str(xi_bloc)];
disp(['<html><p>',str,'</p></html>'])


disp(fprintf('Duration: %s secondes  (code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                                GUI.BLAST_Object.trig_start_bloc));



% Build Figure
% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration*4;
% focus

%% 
% *Courbe noir*  : temps de réaction des réponses en ms;
% *Courbe verte* : niveau local d'instabilité des réponses exprimées en
% pourcentage;
% *Courbe rose*  : niveau local de contamination des anomalies EEG
% critiques ou intercritiques exprimé en pourcentage;
% *Point rouge*  : fausse alarme;
% *Point noir*   : omission;
% *Point bleu*   : bonne réponse;
% *Trait rouge*   : position et durée des anomalies EEG critiques ou intercritiques 


fprintf('\n\n')

%% 
% *Caractéristiques électrophysiologiques pathologiques*
%

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)];

T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other'});
disp(T)
clear datas

fprintf('\n\n\n\n\n\n\n\n')

%% 
% *Cognition : Précision et Vitesse*

% Nb
datas = [round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
    	 round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim];
% pct
datas = [datas, [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss;
                 nan]];
% avg
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw);
                 nan;
                 nan]];
% median
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw);
                 nan;
                 nan]];
             
T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage', 'Avg_ms', 'Med_ms'},...
    'RowNames', {'OK';
                 'False';
                 'Miss';
                 'All'});
disp(T)
clear datas



fprintf('\n\n')
%%
% *Cognition : Pourcentage de Stabilité de l'attention*

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000];


T = array2table(datas,...
    'VariableNames',{'Score'},...
    'RowNames', {'pct20';
                 'pct40';
                 'pct1500';
                 'pct3000'});
disp(T)
clear datas
     
%%
% *PCT 20 et 40*     : stabilité d'attention qui prend en compte la qualité et la
% stabilité des réponses indépendamment de la vitesse des temps de
% réaction; 
% *PCT 1500 et 3000* : stabilité de l'attention qui prend en compte la qualité, la stabilité et la vitesse des réponses.
% 
% *PCT 3000 et 40* : niveau de stabilité de l'attention calculé sur l'ensemble des temps de réaction
% *PCT 1500 et 20* : niveau de stabilité de l'attention calculé sur les temps de réaction les plus rapide uniquement.

fprintf('\n\n')

%%
% *Cognition : Attention focalisée*
% 
datas = [GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max];
     
datas = [datas, [length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)]];
           
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)]];
             
datas = [datas, [round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))]];

             T = array2table(datas,...
    'VariableNames',{'Serie' 'Nb' 'Total_Time', 'Total_Pourc'},...
    'RowNames', {'0-10';
                 '10-20';
                 '20-100'});
disp(T)
clear datas

%%
% *Série*       : Nombre maximum d'item consécutifs qui ont été réalisés à ce niveau
% de performance;
% *Nb*          : Nombre de séries qui ont été réalisées;
% *Total_Time*  : Durée (en seconde) totale que ces séries représentent;
% *Total_Pourc* : Pourcentage du bloc qui a été réalisé à ce niveau de
% performance


fprintf('\n\n\n\n\n\n\n\n\n');





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOC 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% *BLOC 2*

xi_bloc = 2;
if xi_bloc > length(GUI.BLAST_Object.bloc) 
    return
end

str = ['Bloc ', num2str(xi_bloc)];
disp(['<html><p>',str,'</p></html>'])


disp(fprintf('Duration: %s secondes  (code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                                GUI.BLAST_Object.trig_start_bloc));



% Build Figure
% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration*4;
% focus

%% 
% *Courbe noir*  : temps de réaction des réponses en ms;
% *Courbe verte* : niveau local d'instabilité des réponses exprimées en
% pourcentage;
% *Courbe rose*  : niveau local de contamination des anomalies EEG
% critiques ou intercritiques exprimé en pourcentage;
% *Point rouge*  : fausse alarme;
% *Point noir*   : omission;
% *Point bleu*   : bonne réponse;
% *Trait rouge*   : position et durée des anomalies EEG critiques ou intercritiques 


fprintf('\n\n')

%% 
% *Caractéristiques électrophysiologiques pathologiques*
%

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)];

T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other'});
disp(T)
clear datas

fprintf('\n\n\n\n\n\n\n\n')

%% 
% *Cognition : Précision et Vitesse*

% Nb
datas = [round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
    	 round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim];
% pct
datas = [datas, [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss;
                 nan]];
% avg
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw);
                 nan;
                 nan]];
% median
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw);
                 nan;
                 nan]];
             
T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage', 'Avg_ms', 'Med_ms'},...
    'RowNames', {'OK';
                 'False';
                 'Miss';
                 'All'});
disp(T)
clear datas



fprintf('\n\n')
%%
% *Cognition : Pourcentage de Stabilité de l'attention*

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000];


T = array2table(datas,...
    'VariableNames',{'Score'},...
    'RowNames', {'pct20';
                 'pct40';
                 'pct1500';
                 'pct3000'});
disp(T)
clear datas
     
%%
% *PCT 20 et 40*     : stabilité d'attention qui prend en compte la qualité et la
% stabilité des réponses indépendamment de la vitesse des temps de
% réaction; 
% *PCT 1500 et 3000* : stabilité de l'attention qui prend en compte la qualité, la stabilité et la vitesse des réponses.
% 
% *PCT 3000 et 40* : niveau de stabilité de l'attention calculé sur l'ensemble des temps de réaction
% *PCT 1500 et 20* : niveau de stabilité de l'attention calculé sur les temps de réaction les plus rapide uniquement.

fprintf('\n\n')

%%
% *Cognition : Attention focalisée*
% 
datas = [GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max];
     
datas = [datas, [length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)]];
           
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)]];
             
datas = [datas, [round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))]];

             T = array2table(datas,...
    'VariableNames',{'Serie' 'Nb' 'Total_Time', 'Total_Pourc'},...
    'RowNames', {'0-10';
                 '10-20';
                 '20-100'});
disp(T)
clear datas

%%
% *Série*       : Nombre maximum d'item consécutifs qui ont été réalisés à ce niveau
% de performance;
% *Nb*          : Nombre de séries qui ont été réalisées;
% *Total_Time*  : Durée (en seconde) totale que ces séries représentent;
% *Total_Pourc* : Pourcentage du bloc qui a été réalisé à ce niveau de
% performance

fprintf('\n\n\n\n\n\n\n\n\n');





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOC 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% *BLOC 3*

xi_bloc = 3;
if xi_bloc > length(GUI.BLAST_Object.bloc) 
    return
end

str = ['Bloc ', num2str(xi_bloc)];
disp(['<html><p>',str,'</p></html>'])


disp(fprintf('Duration: %s secondes  (code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                                GUI.BLAST_Object.trig_start_bloc));



% Build Figure
% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration*4;
% focus

%% 
% *Courbe noir*  : temps de réaction des réponses en ms;
% *Courbe verte* : niveau local d'instabilité des réponses exprimées en
% pourcentage;
% *Courbe rose*  : niveau local de contamination des anomalies EEG
% critiques ou intercritiques exprimé en pourcentage;
% *Point rouge*  : fausse alarme;
% *Point noir*   : omission;
% *Point bleu*   : bonne réponse;
% *Trait rouge*   : position et durée des anomalies EEG critiques ou intercritiques 


fprintf('\n\n')

%% 
% *Caractéristiques électrophysiologiques pathologiques*
%

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)];

T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other'});
disp(T)
clear datas

fprintf('\n\n\n\n\n\n\n\n')

%% 
% *Cognition : Précision et Vitesse*

% Nb
datas = [round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
    	 round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim];
% pct
datas = [datas, [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss;
                 nan]];
% avg
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw);
                 nan;
                 nan]];
% median
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw);
                 nan;
                 nan]];
             
T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage', 'Avg_ms', 'Med_ms'},...
    'RowNames', {'OK';
                 'False';
                 'Miss';
                 'All'});
disp(T)
clear datas



fprintf('\n\n')
%%
% *Cognition : Pourcentage de Stabilité de l'attention*

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000];


T = array2table(datas,...
    'VariableNames',{'Score'},...
    'RowNames', {'pct20';
                 'pct40';
                 'pct1500';
                 'pct3000'});
disp(T)
clear datas
     
%%
% *PCT 20 et 40*     : stabilité d'attention qui prend en compte la qualité et la
% stabilité des réponses indépendamment de la vitesse des temps de
% réaction; 
% *PCT 1500 et 3000* : stabilité de l'attention qui prend en compte la qualité, la stabilité et la vitesse des réponses.
% 
% *PCT 3000 et 40* : niveau de stabilité de l'attention calculé sur l'ensemble des temps de réaction
% *PCT 1500 et 20* : niveau de stabilité de l'attention calculé sur les temps de réaction les plus rapide uniquement.


fprintf('\n\n')

%%
% *Cognition : Attention focalisée*
% 
datas = [GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max];
     
datas = [datas, [length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)]];
           
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)]];
             
datas = [datas, [round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))]];

             T = array2table(datas,...
    'VariableNames',{'Serie' 'Nb' 'Total_Time', 'Total_Pourc'},...
    'RowNames', {'0-10';
                 '10-20';
                 '20-100'});
disp(T)
clear datas

%%
% *Série*       : Nombre maximum d'item consécutifs qui ont été réalisés à ce niveau
% de performance;
% *Nb*          : Nombre de séries qui ont été réalisées;
% *Total_Time*  : Durée (en seconde) totale que ces séries représentent;
% *Total_Pourc* : Pourcentage du bloc qui a été réalisé à ce niveau de
% performance


fprintf('\n\n\n\n\n\n\n\n\n');





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOC 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% *BLOC 4*

xi_bloc = 4;
if xi_bloc > length(GUI.BLAST_Object.bloc) 
    return
end

str = ['Bloc ', num2str(xi_bloc)];
disp(['<html><p>',str,'</p></html>'])


disp(fprintf('Duration: %s secondes  (code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                                GUI.BLAST_Object.trig_start_bloc));



% Build Figure
% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration*4;
% focus

%% 
% *Courbe noir*  : temps de réaction des réponses en ms;
% *Courbe verte* : niveau local d'instabilité des réponses exprimées en
% pourcentage;
% *Courbe rose*  : niveau local de contamination des anomalies EEG
% critiques ou intercritiques exprimé en pourcentage;
% *Point rouge*  : fausse alarme;
% *Point noir*   : omission;
% *Point bleu*   : bonne réponse;
% *Trait rouge*   : position et durée des anomalies EEG critiques ou intercritiques 


fprintf('\n\n')

%% 
% *Caractéristiques électrophysiologiques pathologiques*
%

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)];

T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other'});
disp(T)
clear datas

fprintf('\n\n\n\n\n\n\n\n')

%% 
% *Cognition : Précision et Vitesse*

% Nb
datas = [round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
    	 round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim];
% pct
datas = [datas, [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss;
                 nan]];
% avg
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw);
                 nan;
                 nan]];
% median
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw);
                 nan;
                 nan]];
             
T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage', 'Avg_ms', 'Med_ms'},...
    'RowNames', {'OK';
                 'False';
                 'Miss';
                 'All'});
disp(T)
clear datas



fprintf('\n\n')
%%
% *Cognition : Pourcentage de Stabilité de l'attention*

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000];


T = array2table(datas,...
    'VariableNames',{'Score'},...
    'RowNames', {'pct20';
                 'pct40';
                 'pct1500';
                 'pct3000'});
disp(T)
clear datas
     
%%
% *PCT 20 et 40*     : stabilité d'attention qui prend en compte la qualité et la
% stabilité des réponses indépendamment de la vitesse des temps de
% réaction; 
% *PCT 1500 et 3000* : stabilité de l'attention qui prend en compte la qualité, la stabilité et la vitesse des réponses.
% 
% *PCT 3000 et 40* : niveau de stabilité de l'attention calculé sur l'ensemble des temps de réaction
% *PCT 1500 et 20* : niveau de stabilité de l'attention calculé sur les temps de réaction les plus rapide uniquement.


fprintf('\n\n')

%%
% *Cognition : Attention focalisée*
% 
datas = [GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max];
     
datas = [datas, [length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)]];
           
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)]];
             
datas = [datas, [round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))]];

             T = array2table(datas,...
    'VariableNames',{'Serie' 'Nb' 'Total_Time', 'Total_Pourc'},...
    'RowNames', {'0-10';
                 '10-20';
                 '20-100'});
disp(T)
clear datas

%%
% *Série*       : Nombre maximum d'item consécutifs qui ont été réalisés à ce niveau
% de performance;
% *Nb*          : Nombre de séries qui ont été réalisées;
% *Total_Time*  : Durée (en seconde) totale que ces séries représentent;
% *Total_Pourc* : Pourcentage du bloc qui a été réalisé à ce niveau de
% performance


fprintf('\n\n\n\n\n\n\n\n\n');





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOC 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% *BLOC 5*

xi_bloc = 5;
if xi_bloc > length(GUI.BLAST_Object.bloc) 
    return
end

str = ['Bloc ', num2str(xi_bloc)];
disp(['<html><p>',str,'</p></html>'])


disp(fprintf('Duration: %s secondes  (code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                                GUI.BLAST_Object.trig_start_bloc));



% Build Figure
% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{xi_bloc}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration*4;
% focus

%% 
% *Courbe noir*  : temps de réaction des réponses en ms;
% *Courbe verte* : niveau local d'instabilité des réponses exprimées en
% pourcentage;
% *Courbe rose*  : niveau local de contamination des anomalies EEG
% critiques ou intercritiques exprimé en pourcentage;
% *Point rouge*  : fausse alarme;
% *Point noir*   : omission;
% *Point bleu*   : bonne réponse;
% *Trait rouge*   : position et durée des anomalies EEG critiques ou intercritiques 


fprintf('\n\n')

%% 
% *Caractéristiques électrophysiologiques pathologiques*
%

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct);
         GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb, round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)];

T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage'},...
    'RowNames', {'AIC Generalized',...
                 'AIC Focal',...
                 'Seizure Generalized',...
                 'Seizure Focal',...
                 'Other'});
disp(T)
clear datas

fprintf('\n\n\n\n\n\n\n\n')

%% 
% *Cognition : Précision et Vitesse*

% Nb
datas = [round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
    	 round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim);...
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim];
% pct
datas = [datas, [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false;...
                 GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss;
                 nan]];
% avg
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw);
                 nan;
                 nan]];
% median
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw);
                 nan;
                 nan]];
             
T = array2table(datas,...
    'VariableNames',{'Nb', 'Pourcentage', 'Avg_ms', 'Med_ms'},...
    'RowNames', {'OK';
                 'False';
                 'Miss';
                 'All'});
disp(T)
clear datas



fprintf('\n\n')
%%
% *Cognition : Pourcentage de Stabilité de l'attention*

datas = [GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500;
         GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000];


T = array2table(datas,...
    'VariableNames',{'Score'},...
    'RowNames', {'pct20';
                 'pct40';
                 'pct1500';
                 'pct3000'});
disp(T)
clear datas
     
%%
% *PCT 20 et 40*     : stabilité d'attention qui prend en compte la qualité et la
% stabilité des réponses indépendamment de la vitesse des temps de
% réaction; 
% *PCT 1500 et 3000* : stabilité de l'attention qui prend en compte la qualité, la stabilité et la vitesse des réponses.
% 
% *PCT 3000 et 40* : niveau de stabilité de l'attention calculé sur l'ensemble des temps de réaction
% *PCT 1500 et 20* : niveau de stabilité de l'attention calculé sur les temps de réaction les plus rapide uniquement.


fprintf('\n\n')

%%
% *Cognition : Attention focalisée*
% 
datas = [GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max;...
         GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max];
     
datas = [datas, [length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur);...
                 length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)]];
           
datas = [datas, [round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur);...
                 round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)]];
             
datas = [datas, [round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration));...
                 round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))]];

             T = array2table(datas,...
    'VariableNames',{'Serie' 'Nb' 'Total_Time', 'Total_Pourc'},...
    'RowNames', {'0-10';
                 '10-20';
                 '20-100'});
disp(T)
clear datas

%%
% *Série*       : Nombre maximum d'item consécutifs qui ont été réalisés à ce niveau
% de performance;
% *Nb*          : Nombre de séries qui ont été réalisées;
% *Total_Time*  : Durée (en seconde) totale que ces séries représentent;
% *Total_Pourc* : Pourcentage du bloc qui a été réalisé à ce niveau de
% performance


fprintf('\n\n\n\n\n\n\n\n\n');




