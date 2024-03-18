function output =  Extract_Stabilo_Vania_Protocol_Bloc_Sample_AIC(TRC_filename, selection)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EX : TRC_filename = '/Volumes/dycog/Epilepto/Vania/Analyse_Clinic/Stabilo/datas_raw/Cognit-AIC_027DR190318.TRC';
%      out =  Extract_Stabilo_Vania_Protocol_Bloc_Sample_AIC(TRC_filename);
%
% MAJ
%   18/11/21    RB
%       add BLAST items selection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global output
output = {};


% Extract trigger
[trig, Fs] = stabilo_read_micromed(TRC_filename);    % en s

% Exctract AIC information
AIC = Extract_AIC_infos(TRC_filename);



ou_start_bloc_95 = find(trig(:,1) == 95);
ou_start_bloc_99 = find(trig(:,1) == 99);


% If there is no 95 code OR LESS THAN 99
% Then we searcgh the 99 code
if (isempty(ou_start_bloc_95) || length(ou_start_bloc_95)<length(ou_start_bloc_99))
    ou_start_bloc = ou_start_bloc_99;
    output.trig_start_bloc = '99';
else
    ou_start_bloc = ou_start_bloc_95;
    output.trig_start_bloc = '95';
end


ou_start_bloc = [ou_start_bloc; length(trig)];

%ou_start_bloc = ou_start_bloc(3:3:end);
% Because some TRC have only two 99 not three 
%ou_start_bloc = ou_start_bloc(find(diff(ou_start_bloc)>1));
% Because some TRC have some 1 between two 95
ou_start_bloc = ou_start_bloc(find(diff(trig(ou_start_bloc,2))>10));
ou_start_bloc = [ou_start_bloc; length(trig)];



ou_stop_bloc = [];
% find last bloc feedback
for xi_start = 1 : length(ou_start_bloc)-1
    
    tmp = find(trig(1:ou_start_bloc(xi_start+1),1) > 70 & trig(1:ou_start_bloc(xi_start+1),1) < 74);
    ou_stop_bloc(xi_start) = tmp(end);
    
end, clear xi_start



output.Fs = Fs;
output.ref_time_session = trig(1,2);   % en s
output.session_duration = trig(end,2) - trig(1,2);
output.session_duration_without_eye = trig(end,2) - trig(ou_start_bloc(1),2);

output.AIC              = AIC;


[output.OY, output.CY] = Extract_Eye(trig);


% find echantilon size minium for all bloc
bloc_size_min = 10000;
bloc_size_max = 0;
for xi_run = 1:length(ou_start_bloc)-1
    code_timing = trig(ou_start_bloc(xi_run):ou_stop_bloc(xi_run), :);
    code_prepare = stabilo_code4JPscore(code_timing);
    if bloc_size_min > size(code_prepare,1)
        bloc_size_min = size(code_prepare,1);
    end
    if size(code_prepare,1) > bloc_size_max  
        bloc_size_max = size(code_prepare,1);
    end
end, clear xi_run


% Selection of Blast items (all or select)
switch selection
    case 'all'
        items_selection.start =1;
        items_selection.nb_items = bloc_size_max;
    case 'select'
        items_selection = UI_specify_items(bloc_size_min, bloc_size_max);
        fprintf('We select %s items, start to %s\n', num2str(items_selection.nb_items), num2str(items_selection.start))
        clear bloc_size_max bloc_size_min
        
end


quel_bloc = 0;
% WARNING il faudra changer le 2
for xi_run = 1:length(ou_start_bloc)-1

    quel_bloc = quel_bloc+1;
    fprintf('Bloc %s\n', num2str(quel_bloc))

    code_timing = trig(ou_start_bloc(xi_run):ou_stop_bloc(xi_run), :);
    
    output.bloc{quel_bloc}.version = Test_Version(code_timing);

    code_prepare = stabilo_code4JPscore(code_timing);

    % Here we select which Blast item we want
    if items_selection.start+items_selection.nb_items-1 > size(code_prepare,1)
        code_prepare = code_prepare([items_selection.start:end],:);
    else
        code_prepare = code_prepare([items_selection.start:items_selection.start+items_selection.nb_items-1],:);
    end
    
    % Infos Blast by stim
    output.bloc{quel_bloc}.sample.Blast.lat = code_prepare(:,1);    % en s
    output.bloc{quel_bloc}.sample.Blast.FB  = code_prepare(:,3);
    output.bloc{quel_bloc}.sample.Blast.RT  = code_prepare(:,4);    % en s
    
    % Extract Blast score for all bloc
    [output.bloc{quel_bloc}.global.Blast,...
     output.bloc{quel_bloc}.sample.Blast.stab] = stabilo_scores_Run_Sampl_AIC([code_prepare(:,5) code_prepare(:,4)]);
    
    start_bloc = min(code_timing(max([code_timing(:,1) == 12, code_timing(:,1) == 11], [], 2), 2));
    stop_bloc  = max(code_timing(max([code_timing(:,1) == 71, code_timing(:,1) == 72, code_timing(:,1) == 73], [], 2), 2));
    output.bloc{quel_bloc}.duration   = stop_bloc - start_bloc;   % de la première stim à la dernière réponse
    clear start_bloc stop_bloc
    output.bloc{quel_bloc}.start_stop = [trig(ou_start_bloc(xi_run),2) trig(ou_stop_bloc(xi_run),2)]; 
        
    % Exctract AIC stat far all bloc
    output.bloc{quel_bloc}.global.AIC = Extract_AIC_all_bloc(AIC.sample(AIC.sample > code_timing(1,2) & AIC.sample < code_timing(end,2)),...
                                                             AIC.label(AIC.sample > code_timing(1,2) & AIC.sample < code_timing(end,2)),...
                                                             quel_bloc);
    
                 
     % AIC / STIM association
     output.bloc{quel_bloc}.sample.AIC = associatio_AIC_STIm(output.bloc{quel_bloc}.sample.Blast);
                     
                     
     output.bloc{quel_bloc}.focus = instability_focus(output.bloc{quel_bloc}.sample.Blast.stab, quel_bloc);
   
                                                         
                                                         
     clear code_timing code_prepare
     
     
end, clear xi_run



%Plot_Bloc_AIC()  



function Plot_Bloc_AIC()  
        
global output

time_aic = output.AIC.sample;
ampl_aic = ones(1,length(time_aic));

time_bloc = [];
ampl_bloc = [];
for xi = 1 : numel(output.bloc)
    time_bloc = [time_bloc output.bloc{xi}.start_stop(1) output.bloc{xi}.start_stop(1) output.bloc{xi}.start_stop(2) output.bloc{xi}.start_stop(2)];
    ampl_bloc = [ampl_bloc 0 2 2 0];
end, clear xi


fill(time_bloc, ampl_bloc, 'k-')
hold on
plot(time_aic, ampl_aic, 'r*')

function sortie = Extract_AIC_infos(filename)

%%%%%%%%%%%%%%%
% Extract AIC label and sample
% from TRC file
%
% sortie : 
%           label  : string, only one letter
%           sample : scalar, latency, seconde
%
%%%%%%%%%%%%%%


sortie.label      = {};
sortie.code_label = [];
sortie.sample     = [];
t = TRC(filename);
for xi = 1 : size(t.m_notes,1)
    
    if strcmpi(deblank(t.m_notes(xi).comment), 'a') || strcmpi(deblank(t.m_notes(xi).comment), 'b') || strcmpi(deblank(t.m_notes(xi).comment), 'c') || strcmpi(deblank(t.m_notes(xi).comment), 'd') ||...
            strcmpi(deblank(t.m_notes(xi).comment), 'e') || strcmpi(deblank(t.m_notes(xi).comment), 'f') || strcmpi(deblank(t.m_notes(xi).comment), 'g') || strcmpi(deblank(t.m_notes(xi).comment), 'h') ||...
            strcmpi(deblank(t.m_notes(xi).comment), 'i') || strcmpi(deblank(t.m_notes(xi).comment), 'j')
        sortie.label{end+1}  = deblank(t.m_notes(xi).comment);
        switch sortie.label{end}
            case 'a'
                sortie.code_label(end+1) = 1;
            case 'b'
                sortie.code_label(end+1) = 1;
            case 'c'
                sortie.code_label(end+1) = 2;
            case 'd'
                sortie.code_label(end+1) = 2;
            case 'e'
                sortie.code_label(end+1) = 3;
            case 'f'
                sortie.code_label(end+1) = 3;
            case 'g'
                sortie.code_label(end+1) = 4;
            case 'h'
                sortie.code_label(end+1) = 4;
            case 'i'
                sortie.code_label(end+1) = 5;
            case 'j'
                sortie.code_label(end+1) = 5;
        end
        
        sortie.sample(end+1) = double(t.m_notes(xi).sample)/double(t.Header.samplingRate);
        
    end
end, clear xi

function  sortie = Extract_AIC_all_bloc(aic_sample, aic_label, quel_bloc)

%%%%%%%%%%%%%%
% 
% AIC 
% a-b : g?n?ralis?e
% c-d : focale
% e-f : autres
% g-h : crise g?n?ralis?, absence
% i-j   : crise focale
%
%%%%%%%%%%%%%

global output


sortie.aic_general_nb   = 0;
sortie.aic_general_pct  = 0;
sortie.aic_focal_nb     = 0;
sortie.aic_focal_pct    = 0;
sortie.other_nb         = 0;
sortie.other_pct        = 0;
sortie.seiz_general_nb  = 0;
sortie.seiz_general_pct = 0;
sortie.seiz_focal_nb    = 0;
sortie.seiz_focal_pct   = 0;


if ~isempty(aic_label)
    % test first instance to know if it's the bengining or the end of the AIC
    switch aic_label{1}
        case 'b'
            start_xi_aic = 2;
            sortie.aic_general_nb = sortie.aic_general_nb + 1;
            sortie.aic_general_pct = sortie.aic_general_pct + ((aic_sample(1) - output.bloc{quel_bloc}.start_stop(1))/output.bloc{quel_bloc}.duration);
        case 'd'
            start_xi_aic = 2;
            sortie.aic_focal_nb = sortie.aic_focal_nb + 1;
            sortie.aic_focal_pct = sortie.aic_focal_pct + ((aic_sample(1) - output.bloc{quel_bloc}.start_stop(1))/output.bloc{quel_bloc}.duration);
        case 'f'
            start_xi_aic = 2;
            sortie.other_nb = sortie.aic_other_nb + 1;
            sortie.other_pct = sortie.other_pct + ((aic_sample(1) - output.bloc{quel_bloc}.start_stop(1))/output.bloc{quel_bloc}.duration);
        case 'h'
            start_xi_aic = 2;
            sortie.seiz_general_nb = sortie.seiz_general_nb + 1;
            sortie.seiz_general_pct = sortie.seiz_general_pct + ((aic_sample(1) - output.bloc{quel_bloc}.start_stop(1))/output.bloc{quel_bloc}.duration);
        case 'j'
            start_xi_aic = 2;
            sortie.seiz_focal_nb = sortie.seiz_focal_nb + 1;
            sortie.seiz_focal_pct = sortie.seiz_focal_pct + ((aic_sample(1) - output.bloc{quel_bloc}.start_stop(1))/output.bloc{quel_bloc}.duration);
        otherwise
            start_xi_aic = 1;
    end
    
    % test last instance to know if it's the bengining or the end of the AIC
    switch aic_label{end}
        case 'a'
            stop_xi_aic = length(aic_label) - 1;
            sortie.aic_general_nb = sortie.aic_general_nb + 1;
            sortie.aic_general_pct = sortie.aic_general_pct + ((output.bloc{quel_bloc}.start_stop(2) - aic_sample(end))/output.bloc{quel_bloc}.duration);
        case 'c'
            stop_xi_aic = length(aic_label) - 1;
            sortie.aic_focal_nb = sortie.aic_focal_nb + 1;
            sortie.aic_focal_pct = sortie.aic_focal_pct + ((output.bloc{quel_bloc}.start_stop(2) - aic_sample(end))/output.bloc{quel_bloc}.duration);
        case 'e'
            stop_xi_aic = length(aic_label) - 1;
            sortie.other_nb = sortie.other_nb + 1;
            sortie.other_pct = sortie.other_pct + ((output.bloc{quel_bloc}.start_stop(2) - aic_sample(end))/output.bloc{quel_bloc}.duration);
        case 'g'
            stop_xi_aic = length(aic_label) - 1;
            sortie.seiz_general_nb = sortie.seiz_general_nb + 1;
            sortie.seiz_general_pct = sortie.seiz_general_pct + ((output.bloc{quel_bloc}.start_stop(2) - aic_sample(end))/output.bloc{quel_bloc}.duration);
        case 'i'
            stop_xi_aic = length(aic_label) - 1;
            sortie.seiz_focal_nb = sortie.seiz_focal_nb + 1;
            sortie.seiz_focal_pct = sortie.seiz_focal_pct + ((output.bloc{quel_bloc}.start_stop(2) - aic_sample(end))/output.bloc{quel_bloc}.duration);
        otherwise
            stop_xi_aic = length(aic_label);
    end
    
    
    for xi_aic = start_xi_aic : 2 : stop_xi_aic
        
        switch aic_label{xi_aic}
            
            case 'a'
                sortie.aic_general_nb = sortie.aic_general_nb + 1;
                sortie.aic_general_pct = sortie.aic_general_pct + ((aic_sample(xi_aic+1) - aic_sample(xi_aic))/output.bloc{quel_bloc}.duration);
            case 'c'
                sortie.aic_focal_nb  = sortie.aic_focal_nb + 1;
                sortie.aic_focal_pct = sortie.aic_focal_pct + ((aic_sample(xi_aic+1) - aic_sample(xi_aic))/output.bloc{quel_bloc}.duration);
            case 'e'
                sortie.other_nb  = sortie.other_nb + 1;
                sortie.other_pct = sortie.other_pct + ((aic_sample(xi_aic+1) - aic_sample(xi_aic))/output.bloc{quel_bloc}.duration);
            case 'g'
                sortie.seiz_general_nb  = sortie.seiz_general_nb + 1;
                sortie.seiz_general_pct = sortie.seiz_general_pct + ((aic_sample(xi_aic+1) - aic_sample(xi_aic))/output.bloc{quel_bloc}.duration);
            case 'i'
                sortie.seiz_focal_nb  = sortie.seiz_focal_nb + 1;
                sortie.seiz_focal_pct = sortie.seiz_focal_pct + ((aic_sample(xi_aic+1) - aic_sample(xi_aic))/output.bloc{quel_bloc}.duration);
        end
    end, clear xi_aic
end

function out_AIC = associatio_AIC_STIm(stim_info)

global output

start_trial_befor_stim = 0.7;  % en s

% all AIC recorded
aic_start = output.AIC.sample(1:2:end-1);     % en s
aic_stop  = output.AIC.sample(2:2:end);       % en s

out_AIC.befor_stim.nb_aic         = nan(1,size(stim_info.FB,1));
out_AIC.befor_stim.cumultime_aic  = nan(1,size(stim_info.FB,1));
out_AIC.befor_stim.cumulPCT_aic   = nan(1,size(stim_info.FB,1));
out_AIC.nb_aic               = zeros(1,size(stim_info.FB,1));
out_AIC.aic_duration         = zeros(1,size(stim_info.FB,1));
out_AIC.aic_pct_duration     = zeros(1,size(stim_info.FB,1));
out_AIC.aic_position_stim    = {};  
out_AIC.type_aic             = {};

for xi_stim = 1 : size(stim_info.RT,1)
    
    % how many aic befor stim ?
    out_AIC.befor_stim.nb_aic(xi_stim)        = sum(aic_start<(stim_info.lat(xi_stim)-start_trial_befor_stim));
    % cumulate time aic befor stim 
    out_AIC.befor_stim.cumultime_aic(xi_stim) = sum(aic_stop(1:out_AIC.befor_stim.nb_aic)' - aic_start(1:out_AIC.befor_stim.nb_aic)');
    % cumulate percent time befor stim
    out_AIC.befor_stim.cumulPCT_aic(xi_stim)  = out_AIC.befor_stim.cumultime_aic(xi_stim)/(stim_info.lat(xi_stim) - output.ref_time_session);
    
    % find aic within stim
    % within if :
    %             - whole the aic within the stim
    %             - just end of the aic within the stim 
    %             - just start of the aic within the stim 
    % stim is the 4 letters display, trial start is 800ms befor 
    quel_tri_aic = find(output.AIC.sample>(stim_info.lat(xi_stim)-start_trial_befor_stim) & output.AIC.sample<(stim_info.lat(xi_stim) + stim_info.RT(xi_stim)));
    
    out_AIC.type_aic{xi_stim} = {};
    % AIC position according to the stim start : -0.8s befor stim (4 letters)
    out_AIC.aic_position_stim{xi_stim} = [];
    is_start = 0;
    
    for xi_aic = 1 : length(quel_tri_aic)
        
        switch output.AIC.label{quel_tri_aic(xi_aic)}
            
            case {'a', 'c', 'e', 'g', 'i'}  
            
                is_start = 1;
                % is AIC is embeded in trial
                if xi_aic < length(quel_tri_aic)                  %& (output.AIC.sample(quel_tri_aic(xi_aic+1)) < (stim_info.lat(xi_stim) + stim_info.RT(xi_stim)))
                    out_AIC.aic_duration(xi_stim) = out_AIC.aic_duration(xi_stim) + output.AIC.sample(quel_tri_aic(xi_aic)+1) - output.AIC.sample(quel_tri_aic(xi_aic));
                else
                    out_AIC.aic_duration(xi_stim) = out_AIC.aic_duration(xi_stim) + (stim_info.lat(xi_stim) + stim_info.RT(xi_stim)) - output.AIC.sample(quel_tri_aic(xi_aic));
                end
                
                out_AIC.nb_aic(xi_stim) = out_AIC.nb_aic(xi_stim) + 1;
                out_AIC.aic_position_stim{xi_stim} = [out_AIC.aic_position_stim{xi_stim}; 
                                                      output.AIC.sample(quel_tri_aic(xi_aic)) - (stim_info.lat(xi_stim)-start_trial_befor_stim),  output.AIC.sample(quel_tri_aic(xi_aic)+1) - (stim_info.lat(xi_stim)-start_trial_befor_stim)];
                                     
                switch output.AIC.label{quel_tri_aic(xi_aic)}
                    case 'a'
                        out_AIC.type_aic{xi_stim}{end+1} = 'aic_general';
                    case 'c'
                        out_AIC.type_aic{xi_stim}{end+1} = 'aic_focal';   
                    case 'e'
                        out_AIC.type_aic{xi_stim}{end+1} = 'other';
                    case 'g'
                        out_AIC.type_aic{xi_stim}{end+1} = 'seiz_general';
                    case 'i'
                        out_AIC.type_aic{xi_stim}{end+1} = 'seiz_focal';
                        
                end
                
            case {'b', 'd', 'f', 'h', 'j'}
                
                % is an AIC start befor ?
                if is_start == 0;
                        out_AIC.aic_duration(xi_stim) = out_AIC.aic_duration(xi_stim) + (output.AIC.sample(quel_tri_aic(xi_aic))-(stim_info.lat(xi_stim)-start_trial_befor_stim));
                        out_AIC.nb_aic(xi_stim) = out_AIC.nb_aic(xi_stim) + 1;
                        out_AIC.aic_position_stim{xi_stim} = [out_AIC.aic_position_stim{xi_stim}; 
                                         output.AIC.sample(quel_tri_aic(xi_aic)-1) - (stim_info.lat(xi_stim)-start_trial_befor_stim),  output.AIC.sample(quel_tri_aic(xi_aic)) - (stim_info.lat(xi_stim)-start_trial_befor_stim)];
               
                switch output.AIC.label{quel_tri_aic(xi_aic)}
                    case 'b'
                        out_AIC.type_aic{xi_stim}{end+1} = 'aic_general';
                    case 'd'
                        out_AIC.type_aic{xi_stim}{end+1} = 'aic_focal';   
                    case 'f'
                        out_AIC.type_aic{xi_stim}{end+1} = 'other';
                    case 'h'
                        out_AIC.type_aic{xi_stim}{end+1} = 'seiz_general';
                    case 'j'
                        out_AIC.type_aic{xi_stim}{end+1} = 'seiz_focal';
                        
                end
                
                end

        end
        
    end
    
    out_AIC.aic_pct_duration(xi_stim) = out_AIC.aic_duration(xi_stim)/((stim_info.lat(xi_stim) + stim_info.RT(xi_stim)) - (stim_info.lat(xi_stim)-start_trial_befor_stim));
    
end

function [OY, CY] = Extract_Eye(trig)


trig_OY = find(trig(:,1) == 121);
OY = {};
OY.seconde = [trig(max(trig_OY),2),...
              trig(max(trig_OY)-1+min(find(trig(max(trig_OY):end,1) == 120)),2)];
OY.sample  = [trig(max(trig_OY),3),...
              trig(max(trig_OY)-1+min(find(trig(max(trig_OY):end,1) == 120)),3)];

  
trig_CY = find(trig(:,1) == 122);
CY = {};
CY.seconde = [trig(max(trig_CY),2),...
              trig(max(trig_CY)-1+min(find(trig(max(trig_CY):end,1) == 120)),2)];
CY.sample =  [trig(max(trig_CY),3),...
              trig(max(trig_CY)-1+min(find(trig(max(trig_CY):end,1) == 120)),3)];

function version = Test_Version(trig)
        
unelettre    = min(find(trig(:,1) == 11));
quatrelettre = unelettre - 1 + min(find(trig(unelettre:end,1) > 100));
ISI = trig(quatrelettre,2) - trig(unelettre,2);

if ISI > 0.8
    version = '-10ans';
else 
    version = '+10ans';
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



