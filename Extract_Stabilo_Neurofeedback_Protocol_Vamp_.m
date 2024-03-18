function output = Extract_Stabilo_Neurofeedback_Protocol_Vamp_(filepath_global)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extract Vamp file from Neurofeedback protocol.
%
% EX:
%       filepath_global = '/Volumes/Backup
%       DD/From_BAIE/Vania/Analyse_Clinic/Stabilo/datas_raw/Vamp/GA299_S01/GA299_blast.ahdr'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global output
output = {};


% LOAD Vamp file
[filepath, filename, ~ ] = fileparts(filepath_global); 
fprintf("Import Vamp file:\n\t\t%s\n", fullfile(filepath, [filename, '_0004.ahdr']))
% Timing in sample
[HDR_vamp, EVT_vamp, ~] = ahdr2mat(fullfile(filepath, [filename, '_0004.ahdr']));
% Move to ms resolution to match with log presentation file
EVT_vamp(:,1) = (EVT_vamp(:,1) / HDR_vamp.Fs)*1000;


% LOAD log presentation file
fprintf("Import log Presentation file:\n\t\t%s\n", fullfile(filepath, [filename, '_S01-mla_mara100.log']))
% Timing in ms
trig = stabilo_read_presentation_log(fullfile(filepath, [filename, '_S01-mla_mara100.log']));
trig(:,2) = trig(:,2)/1000;


% Find timing of the 1st BLAST event in Vamp file
start_Blast_vamp = find(EVT_vamp(:,2) == 15);
start_Blast_vamp = EVT_vamp(start_Blast_vamp(end)+1,1);

% Find timing of the 1st BLAST event in Presentation file
start_Blast_log = find(~isnan(trig(:,2)));
start_Blast_log = trig(start_Blast_log(1), 2);

% 
trig(:,2) = trig(:,2) + (start_Blast_vamp - start_Blast_log);

clear start_Blast_log start_Blast_vamp

    

output.trig_start_bloc = '';
output.Fs = HDR_vamp.Fs;
output.ref_time_session = EVT_vamp(1,1)/output.Fs;   % en s
output.session_duration = HDR_vamp.nSamples/output.Fs;
output.session_duration_without_eye = output.session_duration;

output.AIC.label = {};
output.AIC.code_label = [];
output.AIC.sample = [];
output.code_label = [];
output.OY.seconde = [];
output.OY.sample = [];
output.CY.seconde = [];
output.CY.sample = [];


% dans le protocol Neurofeedback il n'y a qu'un seul bloc
quel_bloc = 1;
code_prepare = stabilo_code4JPscore(trig);

% Selection of Blast items (all or select)
items_selection = UI_specify_items(size(code_prepare,1));
fprintf('We select %s items, start to %s\n', num2str(items_selection.nb_items), num2str(items_selection.start))
clear Bloc_size
code_prepare = code_prepare([items_selection.start:items_selection.start+items_selection.nb_items-1],:);
clear items_selection

output.bloc{quel_bloc}.version = 'Neurofeedback';


% Infos Blast by stim
output.bloc{quel_bloc}.sample.Blast.lat = code_prepare(:,1);    % en s
output.bloc{quel_bloc}.sample.Blast.FB  = code_prepare(:,3);
output.bloc{quel_bloc}.sample.Blast.RT  = code_prepare(:,4);    % en s

% Extract Blast score for all bloc
[output.bloc{quel_bloc}.global.Blast,...
    output.bloc{quel_bloc}.sample.Blast.stab] = stabilo_scores_Run_Sampl_AIC([code_prepare(:,5) code_prepare(:,4)]);

    
code_timing = trig;
    start_bloc = min(code_timing(max([code_timing(:,1) == 12, code_timing(:,1) == 11], [], 2), 2));
    stop_bloc  = max(code_timing(max([code_timing(:,1) == 71, code_timing(:,1) == 72, code_timing(:,1) == 73], [], 2), 2));
    output.bloc{quel_bloc}.duration   = stop_bloc - start_bloc;   % de la première stim à la dernière réponse
    clear start_bloc stop_bloc
    output.bloc{quel_bloc}.start_stop = [min(trig(:,2)) max(trig(:,2))]; 
        
output.bloc{quel_bloc}.global.AIC.aic_general_nb = 0;
output.bloc{quel_bloc}.global.AIC.aic_general_pct = 0;
output.bloc{quel_bloc}.global.AIC.aic_focal_nb = 0;
output.bloc{quel_bloc}.global.AIC.aic_focal_pct = 0;
output.bloc{quel_bloc}.global.AIC.other_nb = 0;
output.bloc{quel_bloc}.global.AIC.other_pct = 0;
output.bloc{quel_bloc}.global.AIC.seiz_general_nb = 0;
output.bloc{quel_bloc}.global.AIC.seiz_general_pct = 0;
output.bloc{quel_bloc}.global.AIC.seiz_focal_nb = 0;
output.bloc{quel_bloc}.global.AIC.seiz_focal_pct = 0;
    
output.bloc{quel_bloc}.sample.AIC.befor_stim.nb_aic = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.befor_stim.cumultime_aic = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.befor_stim.cumulPCT_aic = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.nb_aic = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.aic_duration = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.aic_pct_duration = zeros(1,length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.aic_position_stim = cell(length(code_prepare));
output.bloc{quel_bloc}.sample.AIC.type_aic = cell(length(code_prepare));
    
    
     output.bloc{quel_bloc}.focus = instability_focus(output.bloc{quel_bloc}.sample.Blast.stab, quel_bloc);


     
     
end

function insta_focus_out = instability_focus(instability,quel_bloc)

global output


start_trial_befor_stim = 0.8;  % en s

m10  = int16(instability <= 10);
m20  = int16((instability <= 20) & (instability > 10));
m100 = int16((instability < 100) & (instability > 20));

for xi_stim = 2 : length(m10)-1
    
    if m10(xi_stim) ~= 0
      m10(xi_stim) = m10(xi_stim) + m10(xi_stim-1);
    end
    if m20(xi_stim) ~= 0
      m20(xi_stim) = m20(xi_stim) + m20(xi_stim-1);
    end
    if m100(xi_stim) ~= 0
      m100(xi_stim) = m100(xi_stim) + m100(xi_stim-1);
    end
    
end, clear xi_stim


% how many max successive OK stim 
insta_focus_out.m10.serie_max  = max(m10);
insta_focus_out.m20.serie_max  = max(m20);
insta_focus_out.m100.serie_max = max(m100);

% first successive OK stim 
if max(m10) > 0
    insta_focus_out.m10.lat  = find(m10==max(m10))-(double(max(m10))-1);    % en stim
    % Series duration
    for xi = 1 : length(insta_focus_out.m10.lat)
        insta_focus_out.m10.serie_dur(xi) = (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m10.lat(xi)+(insta_focus_out.m10.serie_max-1)) + output.bloc{quel_bloc}.sample.Blast.RT(insta_focus_out.m10.lat(xi)+(insta_focus_out.m10.serie_max-1))) - (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m10.lat(xi))-start_trial_befor_stim);
    end, clear xi
    insta_focus_out.m10.total_dur = sum(insta_focus_out.m10.serie_dur);
else
    insta_focus_out.m10.lat = [];
    insta_focus_out.m10.total_dur = 0;
    insta_focus_out.m10.serie_dur = [];
end

if max(m20) > 0
    insta_focus_out.m20.lat  = find(m20==max(m20))-(double(max(m20))-1);
    for xi = 1 : length(insta_focus_out.m20.lat)
        insta_focus_out.m20.serie_dur(xi) = (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m20.lat(xi)+(insta_focus_out.m20.serie_max-1)) + output.bloc{quel_bloc}.sample.Blast.RT(insta_focus_out.m20.lat(xi)+(insta_focus_out.m20.serie_max-1))) - (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m20.lat(xi))-start_trial_befor_stim);
    end, clear xi
    insta_focus_out.m20.total_dur = sum(insta_focus_out.m20.serie_dur);
else
    insta_focus_out.m20.lat = [];
    insta_focus_out.m20.total_dur = 0;
    insta_focus_out.m20.serie_dur = [];
end

if max(m100) > 0
    insta_focus_out.m100.lat = find(m100==max(m100))-(double(max(m100))-1);
    for xi = 1 : length(insta_focus_out.m100.lat)
        insta_focus_out.m100.serie_dur(xi) = (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m100.lat(xi)+(insta_focus_out.m100.serie_max-1)) + output.bloc{quel_bloc}.sample.Blast.RT(insta_focus_out.m100.lat(xi)+(insta_focus_out.m100.serie_max-1))) - (output.bloc{quel_bloc}.sample.Blast.lat(insta_focus_out.m100.lat(xi))-start_trial_befor_stim);
    end, clear xi
    insta_focus_out.m100.total_dur = sum(insta_focus_out.m100.serie_dur);
else
    insta_focus_out.m100.lat = [];
    insta_focus_out.m100.total_dur = 0;
    insta_focus_out.m100.serie_dur = [];
end





end