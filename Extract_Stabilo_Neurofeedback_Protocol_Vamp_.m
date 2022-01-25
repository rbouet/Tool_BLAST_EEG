function output = Extract_Stabilo_Neurofeedback_Protocol_Vamp_(filepath_global)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extract Vamp file from Neurofeedback protocol.
%
% EX:
%       filepath_global = '/Volumes/Backup
%       DD/From_BAIE/Vania/Analyse_Clinic/Stabilo/datas_raw/Vamp/GA299_S01/GA299_blast.ahdr'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = {};

[filepath, filename, ~ ] = fileparts(filepath_global); 
[HDR, EVT, DATA] = ahdr2mat(fullfile(filepath, [filename, '_0004.ahdr']));

% a = stabilo_read_presentation_log(fullfile(filepath, [filename, '_S01-mla_mara100.log']))
% 
% % problème, nous n'avons pas le même nombre d'events dans le fichier vamp
% % et dans le log de présentation
% % problème, le timing de présentation débutte avant le timing de vamp
% EVT(1,1)
% a(1,1)
% % Donc les deux fichiers n'ont pas la même référence temporelle
% % Comment les fusionner ?
% 
% % Blast events was wrong, unusual
% % Here the matching:
% % 1 -> 1
% % 2 -> 2
% % 11 -> 11
% % 12 -> 12
% % 101 -> 5
% % 105 -> 9
% % 71 -> 7
% % 72 -> 8
% % 73 -> 9
% EVT(EVT(:,2) == 5, 2) = 101; 
% EVT(EVT(:,2) == 6, 2) = 102; 
% EVT(EVT(:,2) == 7, 2) = 103; 
% EVT(EVT(:,2) == 8, 2) = 104; 
% 
% % on a un problème de doublon ici
% EVT(EVT(:,2) == 9, 2) = 105; 
% EVT(EVT(:,2) == 7, 2) = 71; 
% EVT(EVT(:,2) == 8, 2) = 72; 
% EVT(EVT(:,2) == 9, 2) = 73; 
% 
% 
% EVT = EVT(ismember(EVT(:,2), [1,2,11,12,101:105,71,72,73]),:)
% trig = [EVT(:,2), EVT(:,1)];
% 
% Fs = HDR.Fs;



output.trig_start_bloc = '';
output.Fs = HDR.Fs;
output.ref_time_session = 0;   % en s
output.session_duration = HDR.nSamples/output.Fs;
output.session_duration_without_eye = 0;
output.AIC.label = {};
output.AIC.code_label = [];
output.AIC.sample = [];
output.code_label = [];
output.sample = [];
output.OY.second = [];
output.OY.sample = [];
output.CY.second = [];
output.CY.sample = [];

output.bloc{1}.version = '+10ans';
output.bloc{1}.sample.Blast.lat = [1:5];    % en s
output.bloc{1}.sample.Blast.FB  = [71 71 71 71 71];
output.bloc{1}.sample.Blast.RT  = rand(1,5);    % en s
output.bloc{1}.sample.Blast.stab = [0 0 0 0 0];

output.bloc{1}.global.Blast.nb_stim = 5;
output.bloc{1}.global.Blast.pct_error = 0;
output.bloc{1}.global.Blast.pct_ok = 0;
output.bloc{1}.global.Blast.pct_miss = 0;
output.bloc{1}.global.Blast.pct_false = 0;
output.bloc{1}.global.Blast.avg_raw = 0;
output.bloc{1}.global.Blast.med_raw = 0;
output.bloc{1}.global.Blast.avg_only_good = 0;
output.bloc{1}.global.Blast.med_only_good = 0;
output.bloc{1}.global.Blast.avg_all_plus_penalization = 0;
output.bloc{1}.global.Blast.pct20 = 0;
output.bloc{1}.global.Blast.pct40 = 0;
output.bloc{1}.global.Blast.pct1500 = 0;
output.bloc{1}.global.Blast.pct3000 = 0;

output.bloc{1}.duration = output.session_duration;   % de la première stim à la dernière réponse
output.bloc{1}.start_stop = [0 output.session_duration]; 
output.bloc{1}.global.AIC.aic_general_nb = 0;
output.bloc{1}.global.AIC.aic_general_pct = 0;
output.bloc{1}.global.AIC.aic_focal_nb = 0;
output.bloc{1}.global.AIC.aic_focal_pct = 0;
output.bloc{1}.global.AIC.other_nb = 0;
output.bloc{1}.global.AIC.other_pct = 0;
output.bloc{1}.global.AIC.seiz_general_nb = 0;
output.bloc{1}.global.AIC.seiz_general_pct = 0;
output.bloc{1}.global.AIC.seiz_focal_nb = 0;
output.bloc{1}.global.AIC.seiz_focal_pct = 0;
    
output.bloc{1}.sample.AIC.befor_stim.nb_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.befor_stim.cumultime_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.befor_stim.cumulPCT_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.nb_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_duration = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_pct_duration = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_position_stim = {[],[],[],[],[]};
output.bloc{1}.sample.AIC.type_aic = {[],[],[],[],[]};
     
output.bloc{1}.focus.m10.serie_max = 1;
output.bloc{1}.focus.m10.lat = 1;
output.bloc{1}.focus.m10.serie_dur = 1;
output.bloc{1}.focus.m10.total_dur = 1;
output.bloc{1}.focus.m20.serie_max = 1;
output.bloc{1}.focus.m20.lat = 2;
output.bloc{1}.focus.m20.serie_dur = 1;
output.bloc{1}.focus.m20.total_dur = 1;
output.bloc{1}.focus.m100.serie_max = 1;
output.bloc{1}.focus.m100.lat = 3;
output.bloc{1}.focus.m100.serie_dur = 1;
output.bloc{1}.focus.m100.total_dur = 1;
   
