function GUI_BLAST_Explore()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MAJ
%     08/06/17    copy de header file pour resoudre le probleme de path provenant d'un autre systeme lors de l'import
%     RB
%
%     13/02/18    Modification de l'importation des TRC
%                 Florian sipp class
%     RB
%
%     18/09/18    Modification de l'affichage des marquages
%                 possibilite de fixer le marquage finale
%     RB
%
%     05/03/19    Warning message when no display bandwith was specified
%     RB
%
%     13/05/19    Several channels selection for detection
%     RB
%
%     20/09/20    Add the automatic rapport publication
%     RB
%    
%     16/10/20    Add theta/beta ratio
%     RB
%
%     18/11/21    BLAST's items selection
%     RB
%
%   - 25/01/22   
%               - import neurofeedback Vamp file
%               - create function to create cfg structure for FT structure
%     RB
%
%   - 28/02/22
%               - import BLAST events fom Neurofeedback protocol
%               - Global resum without eyes
%     RB
%
%   - 14/03/24
%               - import multi-files TRC
%                 save .mat object
%                 save .txt file resum all subjects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construction de l'interface 

global GUI
GUI           = {};
GUI.Colors    = [38 94 129;
                 77 144 153;
                 31 57 72;
                 216 212 177]/255;

GUI.Global = {};
GUI.zoom = {};

GUI.Global.fig = figure('Color', GUI.Colors(1,:), 'MenuBar', 'no', 'position',[100 0 1400 850],...
                  'Tag', 'Main_fig', 'CloseRequestFcn',@Main_FigureClose);

Menu_Fil = uimenu('Label', 'File');
           Menu_Clinic = uimenu(Menu_Fil, 'Label', 'Import Clinic');
                         Menu_Clinic_TRC = uimenu(Menu_Clinic, 'Label', 'Explore One TRC', 'Callback', {@Import_TRC});
                         Menu_Clinic_TRC = uimenu(Menu_Clinic, 'Label', 'Extract Many TRC', 'Callback', {@ExtractManyTRC});
           Menu_NeuroFeedback = uimenu(Menu_Fil, 'Label', 'Import NeuroFeedback');
                                uimenu(Menu_NeuroFeedback, 'Label', 'Import matlab', 'Callback', {@Import_Mat_Manu});
                                uimenu(Menu_NeuroFeedback, 'Label', 'Import Vamp', 'Callback', {@Import_Vamp});
           uimenu(Menu_Fil, 'Label', 'Load BLAST_AIC File', 'Callback', {@Load_Blast_Obj});
           uimenu(Menu_Fil, 'Label', 'Save BLAST_AIC file', 'Callback', {@Save_Blast_Obj});

Menu_Dis = uimenu('Label', 'Display');
           uimenu(Menu_Dis, 'Label', 'Global Resum', 'Callback', {@Plot_Bloc_AIC});
           uimenu(Menu_Dis, 'Label', 'Output Text', 'Callback', {@Write_file_output});
Menu_Rap = uimenu('Label', 'Rapport');
           uimenu(Menu_Rap, 'Label', 'Publish Rapport', 'Callback', {@Rapport});
Menu_TF = uimenu('Label', 'TF');
          uimenu(Menu_TF, 'Label', 'Theta/Beta Ratio', 'Callback', {@Launch_RTB});
          uimenu(Menu_TF, 'Label', 'Extract Multi RTB Output', 'Callback', {@Extract_Multi_RTB});
Menu_Modify = uimenu('Label', 'Modify');
              uimenu(Menu_Modify, 'Label', 'Blocs Selection', 'Callback', {@Change_bloc_display});
           
    
           
% Select BLOC 
uicontrol('style','text', 'string', 'Bloc',...
    'position', [30 810 50 30], 'HorizontalAlignment', 'left',...
    'BackgroundColor', GUI.Colors(1,:), 'FontSize', 15);
GUI.Global.Bouton.choix_bloc_plot = uicontrol('Style','list','String', 'None',...
            'max', 50, 'min', 1, 'Tag', 'choix_bloc_plot', ... 
            'position',[20 720 60 100], 'BackgroundColor',GUI.Colors(4,:),...
            'Callback', {@Plot_info_refresh});
           
        
        

GUI.Global.Bouton.Select_plot_info.Grid = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                      'position',[10 680 100 30], 'string', 'Grid', 'tag', 'unselect_grid',...
                                      'FontSize', 12, 'ForegroundColor', [0 0 0],...
                                      'Value', 1, 'Callback', {@Select_display});        
GUI.Global.Bouton.Select_plot_info.RT   = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                      'position',[10 650 100 30], 'string', 'RT', 'tag', 'unselect_RT',...
                                      'FontSize', 12, 'ForegroundColor', [1 1 1],...
                                      'Value', 1, 'Callback', {@Select_display});        
GUI.Global.Bouton.Select_plot_info.FB   = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                      'position',[10 620 100 30], 'string', 'FB', 'tag', 'unselect_FB',...
                                      'FontSize', 12, 'ForegroundColor', [0 0 0],...
                                      'Value', 1, 'Callback', {@Select_display});        
GUI.Global.Bouton.Select_plot_info.Unstability = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                             'position',[10 590 100 30], 'string', 'Unstability', 'tag', 'unselect_unstability',...
                                             'FontSize', 12, 'ForegroundColor', [.35 0.7 0],...
                                             'Value', 1, 'Callback', {@Select_display});        
GUI.Global.Bouton.Select_plot_info.Pct_contamin = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                              'position',[10 560 120 30], 'string', 'Contamination', 'tag', 'unselect_contamin',...
                                              'FontSize', 12, 'ForegroundColor', [.83 0 1],...
                                              'Value', 1, 'Callback', {@Select_display});        
GUI.Global.Bouton.Select_plot_info.focus = uicontrol('Style', 'radiobutton','BackgroundColor', GUI.Colors(1,:),...
                                              'position',[10 530 100 30], 'string', 'Focus', 'tag', 'unselect_focus',...
                                              'FontSize', 12, 'ForegroundColor', [0 0 0],...
                                              'Value', 1, 'Callback', {@Select_display});        

                                          
% Zoom
GUI.Global.Bouton.Zoom = uicontrol('Style', 'pushbutton', 'String', 'Zoom',...
                            'Position', [10 100 100 30], 'BackgroundColor', GUI.Colors(2,:), 'tag', 'Zoome',...
                            'FontSize', 16, 'Callback', {@Zoom_plot}); 

% ClicClac
GUI.Global.Bouton.ClicClac = uicontrol('Style', 'pushbutton', 'String', 'ClicClac',...
                            'Position', [1100 810 100 30], 'BackgroundColor', GUI.Colors(2,:), 'tag', 'ClicClac_Global',...
                            'FontSize', 16, 'Callback', {@ClicClac}); 

                        
                       
 
% normalisation du GUI
child = GUI.Global.fig.Children;
for xi = 1 : length(child)
    if ~strcmp(child(xi).Type, 'uimenu')
        set(child(xi), 'Units', 'normalized');
    end
end

Screen_size = get(0,'ScreenSize');
GUI.Global.fig.Position = [0 0 Screen_size(3)-10 Screen_size(4)-10];



function Main_FigureClose(~,~)
%close all figures
str_CloseSel = questdlg('Are you sure you want to close all windows?',...
    'Close Request Dialog','Yes','No','Yes');

switch str_CloseSel,
    case 'Yes',
        
        delete(findobj('Tag', 'Main_fig'))
        clc;clear;close all; clear global
        
    case 'No'
        return
end
            
function Import_TRC(hObj,evnt)
    
global GUI
[file_name file_path file_ext] = uigetfile({'*.TRC', 'Select TRC file'},...
                                            'MultiSelect', 'off');
                           
fprintf('Extract BLAST scores and AIC...\n')
GUI.BLAST_Object =  Extract_Stabilo_Vania_Protocol_Bloc_Sample_AIC(fullfile(file_path, file_name), 'select');
GUI.BLAST_Object.Path_TRC = [file_path file_name];
GUI.Source = 'TRC_Blast';

% ajust liste des blocs
ob.cb = findobj('Tag', 'choix_bloc_plot');
ob.cb.String = strsplit(num2str([1:size(GUI.BLAST_Object.bloc,2)]));

Plot_info_first()

function Import_Vamp(hObj,evnt)

global GUI
[file_name file_path file_ext] = uigetfile({'*.ahdr', 'Select V-amp file'},...
                                            'MultiSelect', 'off');
                                        
% There are 3 reccords (blast, oe and ce)
% here we specify the gneric filename
file_name = file_name(1:strfind(file_name, '_')-1);


fprintf('There is only 1  BLAST''s bloc and NO AIC.\n')
% en construction, il manque les code correspondant au blast
GUI.BLAST_Object =  Extract_Stabilo_Neurofeedback_Protocol_Vamp_(fullfile(file_path, file_name));
GUI.BLAST_Object.Path_TRC = fullfile(file_path, file_name);
GUI.Source = 'Vamp_Neurofeedback';

% ajust liste des blocs
ob.cb = findobj('Tag', 'choix_bloc_plot');
ob.cb.String = '1';

Plot_info_first()

function Import_Mat_Manu(hObj,evnt)
    
global GUI
[file_name file_path file_ext] = uigetfile({'*.mat', 'Select MAT Manu file'},...
                                            'MultiSelect', 'off');

fprintf('There is no BLAST scores and AIC.\n')
GUI.BLAST_Object =  Extract_Stabilo_Manu_Protocol_Bloc_Sample_AIC(fullfile(file_path, file_name));
GUI.BLAST_Object.Path_TRC = [file_path file_name];
GUI.Source = 'MAT_Neurofeedback';

% ajust liste des blocs
ob.cb = findobj('Tag', 'choix_bloc_plot');
ob.cb.String = '1';

Plot_info_first()

function Save_Blast_Obj(hObj,evnt)
    
global GUI

BLAST_Object = GUI.BLAST_Object;

[TRC_pathname, TRC_filename] = fileparts(GUI.BLAST_Object.Path_TRC);
[filename, pathname] = uiputfile([fullfile(TRC_pathname, TRC_filename) '_BLAST_Obj.mat'], 'Save BLAST Object as');

save(fullfile(pathname, filename), 'BLAST_Object')

function Load_Blast_Obj(hObj,evnt)
    
global GUI

[file_name file_path file_ext] = uigetfile({'*.mat', 'Select BLAST_Object file'},...
                                            'MultiSelect', 'off');
load(fullfile(file_path, file_name))

GUI.BLAST_Object = BLAST_Object;

% ajust liste des blocs
ob.cb = findobj('Tag', 'choix_bloc_plot');
%ob.cb.String = num2strcell([1:size(GUI.BLAST_Object.bloc,2)]);
ob.cb.String = strsplit(num2str([1:size(GUI.BLAST_Object.bloc,2)]))

Plot_info_first()


function ExtractManyTRC(hObj,evnt)
    
global GUI
[file_name file_path file_ext] = uigetfile({'*.TRC', 'Select TRC file'},...
                                            'MultiSelect', 'on');

fid = fopen(fullfile(file_path , ['Resume_', date, '_', num2str(numel(file_name)), 'TRC.txt']), 'w+');
fprintf('save file:\t\t%s\n', fullfile(file_path , ['Resume_', date, '_', num2str(numel(file_name)), 'TRC.txt']))
fprintf(fid, 'Subject\tVersion\tBloc\tStim\t');
fprintf(fid, 'Blast_Bloc_Pct_Error\tBlast_Bloc_Pct_ok\tBlast_Bloc_Pct_Miss\tBlast_Bloc_Pct_False\tBlast_Bloc_TR_Good_Mean\tBlast_Bloc_TR_All_Mean\tBlast_Bloc_PCT20\tBlast_Bloc_PCT40\tBlast_Bloc_PCT1500\tBlast_Bloc_PCT300\t');
fprintf(fid, 'AIC_Bloc_Nb\tAIC_Bloc_Pct\t');
fprintf(fid, 'AIC_Stim_Nb\tAIC_Stim_Length\tAIC_Stim_Pct\tAIC_Stim_Befor_Nb\tAIC_Stim_Befor_Length_Cumul\t');
fprintf(fid, 'Blast_Stim_FB\tBlast_Stim_RT\n');

% for all TRC files
for xi_trc = 1:numel(file_name)
    % TRC Extraction
    fprintf('Extraction file:\t\t%s\n', file_name{xi_trc})
    BLAST_Object =  Extract_Stabilo_Vania_Protocol_Bloc_Sample_AIC(fullfile(file_path, file_name{xi_trc}), 'all');
    BLAST_Object.Path_TRC = [file_path file_name{xi_trc}];

    % Save table
    for xi_bloc = 1 : numel(BLAST_Object.bloc)
        for xi_stim= 1 : length(BLAST_Object.bloc{xi_bloc}.sample.AIC.nb_aic)
            fprintf(fid, '%s\t%s\t%s\t', file_name{xi_trc}, BLAST_Object.bloc{xi_bloc}.version, num2str(xi_bloc), num2str(xi_stim));
            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct_error),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.avg_all_plus_penalization),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct20),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct40),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500),...
                num2str(BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000));
            fprintf(fid, '%s\t%s\t',...
                num2str(sum(BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb + BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb +...
                BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb + BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_nb +...
                BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb)),...
                num2str(sum(BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct + BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct +...
                BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct + BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct +...
                BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct)));
            
            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t',...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.AIC.nb_aic(xi_stim)),...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_duration(xi_stim)),...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.AIC.aic_pct_duration(xi_stim)),...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.AIC.befor_stim.nb_aic(xi_stim)),...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.AIC.befor_stim.cumultime_aic(xi_stim)));
            fprintf(fid, '%s\t%s\n',...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.Blast.FB(xi_stim)),...
                num2str(BLAST_Object.bloc{xi_bloc}.sample.Blast.RT(xi_stim)));
            
        end % xi_stim
        clear xi_stim
        
    end % xi_bloc    
    clear xi_bloc
    
    % Save .mat
    [~, filename, ~] = fileparts(file_name{xi_trc});
    filename = [filename, '_BLAST_Object', '.mat'];
    fprintf('save file:\t\t%s\n', filename)
    save(fullfile(file_path, filename), 'BLAST_Object')
    clear BLAST_Object filename
    
end % xi_trc
fclose(fid);


function Plot_info_first(hObj,evnt)

global GUI

% quel Bloc ?
ob.cb = findobj('Tag', 'choix_bloc_plot');



% ADD another axes to the right
firstAxes = axes('Color', 'none', 'FontSize', 20,...
                'YLim', [-1 4],...
                'YAxisLocation', 'right',...
                'XTick', [],...
                'YTick', [-1:4], 'YTickLabel', {''; '0'; '25'; '50'; '75'; '100'},...
                'Box', 'off',...
                'Tag', 'First_axes');
firstAxes.YLabel.String = 'Percentage %';

hold on
NewAxes = axes('Position', firstAxes.Position,... 
                'Color', 'none', 'FontSize', 20,...
                'YLim', [-1 4],...
                'YAxisLocation', 'left',...
                'XTick', [],...
                'YTick', [-1:4], 'YTickLabel', {'-1'; '0'; '1'; '2'; '3'; '4'},...
                'Box', 'off', 'YGrid', 'on');
NewAxes.YLabel.String = 'Reaction Time (s)';
NewAxes.XLabel.String = 'Time (s)';

hold on

% Perturbation                                   
GUI.Global.fig_Features.Plot_pct_perturb = plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat, GUI.BLAST_Object.bloc{ob.cb.Value}.sample.AIC.aic_pct_duration*100/25   ,...
                                     'color', [.83 0 1], 'LineWidth', 3,  'Tag', 'plot_pct_perturb',...
                                     'Visible', 'on');
% Stability
GUI.Global.fig_Features.Plot_stability =  plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat, GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.stab/25,...
                                    'color', [.35 0.7 0], 'LineWidth', 3, 'Tag', 'plot_stab',...
                                    'Visible', 'on');

  
% ADD focus regions
Add_Focus()

% Plot Rt
GUI.Global.fig_Features.Plot_RT = plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat,...
                                GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
                               'color', [1 1 1], 'LineWidth', 3, ...
                               'Tag', 'plot_RT');

GUI.Global.fig_Features.Plot_RT.Parent.XLabel.FontSize = 20;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(end)+5];
GUI.Global.fig_Features.Plot_RT.Parent.FontSize = 20;
GUI.Global.fig_Features.Plot_RT.Parent.YTickLabel{1} = '';



hold on
% FB
col = zeros(length(GUI.BLAST_Object.bloc{1}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1);   
GUI.Global.fig_Features.Plot_FB = scatter(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat, GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
                                   [], col, 'filled', 'MarkerEdgeColor', 'k', ...
                                   'Visible', 'on', 'SizeData', 100, 'Tag', 'plot_FB');
clear col


% plot AIC 
color_aic = [1 0 0;
             0 1 0;
             0 0 1;
             1 .2 .2;
             .2 .2 1];   
for xi_aic = 1 : 2 : length(GUI.BLAST_Object.AIC.sample)
      GUI.Global.fig_Features.Plot_aic_pos(xi_aic) = fill(GUI.BLAST_Object.AIC.sample([xi_aic xi_aic+1 xi_aic+1 xi_aic]),...
                                                   [-.2 -.2 -.3 -.3],...
                                                    ' ', 'Tag', 'plot_aic_pos');  
      GUI.Global.fig_Features.Plot_aic_pos(xi_aic).FaceColor = color_aic(GUI.BLAST_Object.AIC.code_label(xi_aic),:);
      GUI.Global.fig_Features.Plot_aic_pos(xi_aic).EdgeColor = [0 0 0];
end, clear xi_aic

function Add_Focus()

global GUI

% GUI.Global.fig_Features.Plot_focus_10 = [];
% GUI.Global.fig_Features.Plot_focus_20 = [];
% GUI.Global.fig_Features.Plot_focus_100 = [];

% quel Bloc ?
ob.cb = findobj('Tag', 'choix_bloc_plot');

if GUI.Global.Bouton.Select_plot_info.focus.Value == 1
    Visi = 'on';
else 
    Visi = 'off';
end 


for xi_100 = 1 : length(GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m100.lat)
    x_lat = [GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m100.lat(xi_100):GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m100.lat(xi_100)+GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m100.serie_max-1];
    y_lat = [min(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))-.05 max(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))+.05];
    x_lat = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(x_lat);
    x_lat = [x_lat(1)-1 x_lat(end)+1];
    GUI.Global.fig_Features.Plot_focus_100(xi_100) = fill([x_lat x_lat(end:-1:1)], [y_lat(1) y_lat(1) y_lat(2) y_lat(2)],...
                                                    'r', 'FaceAlpha', .3, 'EdgeAlpha', 0, 'Tag', 'focus_100',...
                                                    'Visible', Visi);
end, clear xi_100 x_lat y_lat

for xi_20 = 1 : length(GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m20.lat)
    x_lat = [GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m20.lat(xi_20):GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m20.lat(xi_20)+GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m20.serie_max-1];
    y_lat = [min(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))-.05 max(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))+.05];
    x_lat = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(x_lat);
    x_lat = [x_lat(1)-1 x_lat(end)+1];
    GUI.Global.fig_Features.Plot_focus_20(xi_20) = fill([x_lat x_lat(end:-1:1)], [y_lat(1) y_lat(1) y_lat(2) y_lat(2)],...
                                                   'y', 'FaceAlpha', .3, 'EdgeAlpha', 0, 'Tag', 'focus_20',...
                                                   'Visible', Visi);
end, clear xi_20 x_lat y_lat

for xi_10 = 1 : length(GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m10.lat)
    x_lat = [GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m10.lat(xi_10):GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m10.lat(xi_10)+GUI.BLAST_Object.bloc{ob.cb.Value}.focus.m10.serie_max-1];
    y_lat = [min(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))-.05 max(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT(x_lat))+.05];
    x_lat = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(x_lat);
    x_lat = [x_lat(1)-1 x_lat(end)+1];
    GUI.Global.fig_Features.Plot_focus_10(xi_10) = fill([x_lat x_lat(end:-1:1)], [y_lat(1) y_lat(1) y_lat(2) y_lat(2)],...
                                                   'g', 'FaceAlpha', .3, 'EdgeAlpha', 0, 'Tag', 'focus_10',...
                                                   'Visible', Visi);
end, clear xi_10 x_lat y_lat

function Del_Focus()

global GUI

if isfield(GUI.Global.fig_Features, 'Plot_focus_100')
    if ~isempty(GUI.Global.fig_Features.Plot_focus_100)
        for xi_100 = 1 : length(GUI.Global.fig_Features.Plot_focus_100)
            delete(GUI.Global.fig_Features.Plot_focus_100(xi_100))
        end, clear xi_100
        GUI.Global.fig_Features = rmfield(GUI.Global.fig_Features, 'Plot_focus_100');
    end
end

if isfield(GUI.Global.fig_Features, 'Plot_focus_20')
    if ~isempty(GUI.Global.fig_Features.Plot_focus_20)
        for xi_20 = 1 : length(GUI.Global.fig_Features.Plot_focus_20)
            delete(GUI.Global.fig_Features.Plot_focus_20(xi_20))
        end, clear xi_20
        GUI.Global.fig_Features = rmfield(GUI.Global.fig_Features, 'Plot_focus_20');
    end
end

if isfield(GUI.Global.fig_Features, 'Plot_focus_10')
    if ~isempty(GUI.Global.fig_Features.Plot_focus_10)
        for xi_10 = 1 : length(GUI.Global.fig_Features.Plot_focus_10)
            delete(GUI.Global.fig_Features.Plot_focus_10(xi_10))
        end, clear xi_10
        GUI.Global.fig_Features = rmfield(GUI.Global.fig_Features, 'Plot_focus_10');
    end
end

function Plot_info_refresh(hObj,evnt)

global GUI

% quel Bloc ?
ob.cb = findobj('Tag', 'choix_bloc_plot');

% RT
GUI.Global.fig_Features.Plot_RT.XData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_RT.YData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT;
GUI.Global.fig_Features.Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(end)+5];
% Feedback
col = zeros(length(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB),3);
col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1);
col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1); 
GUI.Global.fig_Features.Plot_FB.CData = col;
clear col
GUI.Global.fig_Features.Plot_FB.XData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_FB.YData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT;
% Stability
GUI.Global.fig_Features.Plot_stability.XData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_stability.YData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.stab/25;
% Perturbation
GUI.Global.fig_Features.Plot_pct_perturb.XData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat;
GUI.Global.fig_Features.Plot_pct_perturb.YData = GUI.BLAST_Object.bloc{ob.cb.Value}.sample.AIC.aic_pct_duration*4;

% focus
Del_Focus()
Add_Focus()

function Plot_Bloc_AIC(hObj,evnt)
        
global GUI


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

[~, tit, ~] = fileparts(GUI.BLAST_Object.Path_TRC);
title(tit, 'interpreter', 'none')

% ClicClac
GUI.fig_AIC_Pos.ClicClac = uicontrol('Style', 'pushbutton', 'String', 'ClicClac',...
                            'Position', [900 470 100 30], 'BackgroundColor', GUI.Colors(2,:), 'tag', 'ClicClac_Blocs',...
                            'FontSize', 16, 'Callback', {@ClicClac}); 
                        
                        
% normalisation du GUI
child = GUI.fig_AIC_Pos.main.Children;
for xi = 1 : length(child)
    if ~strcmp(child(xi).Type, 'uimenu')
        set(child(xi), 'Units', 'normalized');
    end
end
                        
function Select_display(hObj,evnt)

global GUI

if hObj.Value == 1
    Visi = 'on';
else 
    Visi = 'off';
end 

switch hObj.String
    case 'Grid'
        GUI.Global.fig_Features.Plot_RT.Parent.YGrid = Visi;
    case 'RT'
        GUI.Global.fig_Features.Plot_RT.Visible = Visi;
    case 'FB'
        GUI.Global.fig_Features.Plot_FB.Visible = Visi;
    case 'Unstability'
        GUI.Global.fig_Features.Plot_stability.Visible = Visi;
    case 'Contamination'
        GUI.Global.fig_Features.Plot_pct_perturb.Visible = Visi;
    case 'Focus'
        if isfield(GUI.Global.fig_Features, 'Plot_focus_100')
            for xi_100 = 1 : length(GUI.Global.fig_Features.Plot_focus_100)
                GUI.Global.fig_Features.Plot_focus_100(xi_100).Visible = Visi;
            end, clear xi_100
        end
        if isfield(GUI.Global.fig_Features, 'Plot_focus_20')
            for xi_20 = 1 : length(GUI.Global.fig_Features.Plot_focus_20)
                GUI.Global.fig_Features.Plot_focus_20(xi_20).Visible = Visi;
            end, clear xi_20
        end
        if isfield(GUI.Global.fig_Features, 'Plot_focus_10')
            for xi_10 = 1 : length(GUI.Global.fig_Features.Plot_focus_10)
                GUI.Global.fig_Features.Plot_focus_10(xi_10).Visible = Visi;
            end, clear xi_10
        end

    
    
end

function Zoom_plot(hObj,evnt)

global GUI

[gx, ~] = ginput(2);

% if there is another second
if isempty(GUI.zoom)
    quel_second = 1;
else
    quel_second = length(GUI.zoom) + 1;
end


GUI.zoom(quel_second).main = figure('Color', GUI.Colors(1,:), 'MenuBar', 'no', 'position',[00 50 750 1000],...
                                      'Tag', ['Second_' num2str(quel_second)], 'CloseRequestFcn',@Second_FigureClose,...
                                      'WindowButtonMotionFcn', @MouseCallback);
              
              % quel Bloc ?
ob.cb = findobj('Tag', 'choix_bloc_plot');

% ADD another axes to the right
firstAxes = axes('Color', 'none', 'FontSize', 20,...
                'YLim', [-1 5],...
                'YAxisLocation', 'right',...
                'XTick', [],...
                'YTick', [-1:5], 'YTickLabel', {''; '0'; '20'; '40'; '60'; '80'; '100'},...
                'Box', 'off');
firstAxes.YLabel.String = 'Percentage %';

hold on
NewAxes = axes('Position', firstAxes.Position,... 
                'Color', 'none', 'FontSize', 20,...
                'YLim', [-1 3],...
                'YAxisLocation', 'left',...
                'YTick', [-1:3], 'YTickLabel', {'-1'; '0'; '1'; '2'; '3'},...
                'Box', 'off', 'YGrid', 'on');
NewAxes.YLabel.String = 'Reaction Time (s)';
NewAxes.XLabel.String = 'Time (s)';

hold on

% Plot Rt
GUI.zoom(quel_second).Plot_RT = plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat+GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
                                       GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
                                       'color', [1 1 1], 'LineWidth', 3, ...
                                       'Tag', ['seconde_plot_RT_' num2str(quel_second)]);
GUI.zoom(quel_second).Plot_RT.Parent.XLabel.FontSize = 20;
GUI.zoom(quel_second).Plot_RT.Parent.XLim = [GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(1)-5 GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(end)+5];
GUI.zoom(quel_second).Plot_RT.Parent.FontSize = 20;
GUI.zoom(quel_second).Plot_RT.Parent.YTickLabel{1} = '';

% paradigme
lesY = ones(length(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat),1);
latence_between_target_stim = 0.7;  % en s
plot([GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat]',...
      [-lesY lesY-1]', 'k-',...
      'LineWidth', 2);
plot([GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat]'-latence_between_target_stim,...
      [-lesY lesY-1.5]', 'k:',...
      'LineWidth', 2);
clear lesY latence_between_target_stim
  

hold on
% FB
if GUI.Global.Bouton.Select_plot_info.FB.Value
    
    col = zeros(length(GUI.BLAST_Object.bloc{1}.sample.Blast.FB),3);
    col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1:3) = repmat([0.2863    0.9412    0.7255], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 71),1);
    col(find(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1:3) = repmat([0.8745    0.2392    0.0431], sum(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.FB == 72),1);  
    scatter(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat+GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
            GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.RT,...
            [], col, 'filled', 'MarkerEdgeColor', 'k',...
            'Visible', 'on', 'SizeData', 100);
    clear col
end

% inStability
if GUI.Global.Bouton.Select_plot_info.Unstability.Value
    
    plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat, GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.stab/33.33,...
        'color', [.35 0.7 0], 'LineWidth', 3,...
        'Visible', 'on');
end

% Perturbation                                   
if GUI.Global.Bouton.Select_plot_info.Pct_contamin.Value
    
    plot(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat, GUI.BLAST_Object.bloc{ob.cb.Value}.sample.AIC.aic_pct_duration/.33,...
                                     'color', [.83 0 1], 'LineWidth', 3,...
                                     'Visible', 'on');
end

% plot AIC 
color_aic = [1 0 0;
             0 1 0;
             0 0 1;
             1 .2 .2;
             .2 .2 1];   
for xi_aic = 1 : 2 : length(GUI.BLAST_Object.AIC.sample)
      Plot_aic_pos(xi_aic) = fill(GUI.BLAST_Object.AIC.sample([xi_aic xi_aic+1 xi_aic+1 xi_aic]),...
                                  [-.2 -.2 -.5 -.5], ' ');  
      Plot_aic_pos(xi_aic).FaceColor = color_aic(GUI.BLAST_Object.AIC.code_label(xi_aic),:);
      Plot_aic_pos(xi_aic).EdgeColor = [0 0 0];
     
end, clear xi_aic

GUI.zoom(quel_second).Plot_RT.Parent.XLim = gx;


% Slider
GUI.zoom(quel_second).Slider_time = uicontrol('Style', 'slider', 'Min', GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(1), 'Max', GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat(end), 'Tag', ['slider_time_' num2str(quel_second)],...
                                                'SliderStep', [1/length(GUI.BLAST_Object.bloc{ob.cb.Value}.sample.Blast.lat) 0.05],...
                                                'Value', mean(gx),...
                                                'BackgroundColor',GUI.Colors(3,:), 'Position', [150 20 450 15],...
                                                'Callback', {@Slider_time});
% ADD focus regions
%Add_Focus()

% ClicClac
GUI.zoom(quel_second).Bouton.ClicClac = uicontrol('Style', 'pushbutton', 'String', 'ClicClac',...
                            'Position', [660 960 80 30], 'BackgroundColor', GUI.Colors(2,:), 'tag', ['ClicClac_Zoom_' num2str(quel_second)],...
                            'FontSize', 16, 'Callback', {@ClicClac}); 


% normalisation du GUI
child = GUI.zoom(quel_second).main.Children;
for xi = 1 : length(child)
    if ~strcmp(child(xi).Type, 'uimenu')
        set(child(xi), 'Units', 'normalized');
    end
end

Screen_size = get(0,'ScreenSize');
GUI.zoom(quel_second).Position = [0 0 Screen_size(3)/2 Screen_size(4)-10];

function Second_FigureClose(~,~)

delete(gcf)

function Slider_time(hObj,evnt)

global GUI

ob_plt = findobj('Tag', ['seconde_plot_RT_' hObj.Tag(13:end)]);
ob_plt.Parent.XLim = [hObj.Value - [ob_plt.Parent.XLim(2) - ob_plt.Parent.XLim(1)]/2 hObj.Value + [ob_plt.Parent.XLim(2) - ob_plt.Parent.XLim(1)]/2];

function out = MouseCallback(hObject,obj,eventdata,handles)

quel_axis = gca;
coo = get(quel_axis,'Currentpoint');
out = num2str(coo(1,2));
title(gca, ['RT = ', out]);

function Write_file_output(hObj,evnt)

global GUI

[pathTRC, nameTRC, extTRC] = fileparts(GUI.BLAST_Object.Path_TRC);
[filename, pathname] = uiputfile([fullfile(pathTRC, nameTRC) '_BLAST_Rapport.txt'], 'Save BLAST Rapport as');

fid = fopen(fullfile(pathname, filename), 'w+');
clear filename 



fprintf(fid,'File name\n%s\n\n', fullfile(pathTRC, [nameTRC, extTRC]));
clear pathTRC nameTRC

fprintf(fid, '\t\tEpileptic Events\n');
fprintf(fid, 'AIC Generalized\t\t%s\n', num2str(sum(strcmpi(GUI.BLAST_Object.AIC.label, 'a'))));
fprintf(fid, 'AIC Focal\t\t%s\n', num2str(sum(strcmpi(GUI.BLAST_Object.AIC.label, 'c'))));
fprintf(fid, 'Seizure Generalized\t%s\n', num2str(sum(strcmpi(GUI.BLAST_Object.AIC.label, 'g'))));
fprintf(fid, 'Seizure Focal\t\t%s\n', num2str(sum(strcmpi(GUI.BLAST_Object.AIC.label, 'i'))));
fprintf(fid, 'Other\t\t\t%s\n', num2str(sum(strcmpi(GUI.BLAST_Object.AIC.label, 'e'))));
fprintf(fid, 'Total\t\t\t%s\n', num2str(length(GUI.BLAST_Object.AIC.code_label)/2));
fprintf(fid, '\n\n\n');

fprintf(fid, '\t\tTotal\t\tWithout Eyes blocs\n');
fprintf(fid, 'Duration\t%.1f min\t%.1f min\n', GUI.BLAST_Object.session_duration/60, GUI.BLAST_Object.session_duration_without_eye/60);
fprintf(fid, '\n\n\n');
 
for xi_bloc = 1 : length(GUI.BLAST_Object.bloc)
    fprintf(fid, '______\nBLOC %s\n\n', num2str(xi_bloc));
    fprintf(fid, 'Duration: %s secondes\t\t(code %s)\n\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.duration)),...
                                              GUI.BLAST_Object.trig_start_bloc);
    fprintf(fid, 'Nb AIC Generalized: %s (%s%%)\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_nb),...
                                                    num2str(round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_general_pct)));
    fprintf(fid, 'Nb AIC Focal: %s (%s%%)\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_nb),...
                                              num2str(round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.aic_focal_pct)));
    fprintf(fid, 'Seizure Generalized: %s (%s%%)\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_nb),...
                                                     num2str(round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_general_pct)));
    fprintf(fid, 'Seizure Focal: %s (%s%%)\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct),...
                                               num2str(round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.seiz_focal_pct)));
    fprintf(fid, 'Other: %s (%s%%)\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_nb),...
                                       num2str(round(100*GUI.BLAST_Object.bloc{xi_bloc}.global.AIC.other_pct)));
    fprintf(fid, '\n\t\tOK\tFALSE\tMISS\tALL\n');
    % Nb 
    fprintf(fid, 'Nb\t\t%s\t%s\t%s\t%s\n', num2str(round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim)),...
                                       num2str(round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim)),...
                                       num2str(round((GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss/100)*GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim)),...
                                       num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.nb_stim));
    % Percentage
    fprintf(fid, '%%\t\t%s\t%s\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_ok),...
                                           num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_false),...
                                           num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct_miss));
     % avg
     fprintf(fid, 'Avg (ms)\t%s\t\t\t%s\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_only_good)),...
                                      num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.avg_raw)));                                          
     % median
     fprintf(fid, 'Med (ms)\t%s\t\t\t%s\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_only_good)),...
                                      num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.med_raw)));
    % PCT
    fprintf(fid, '\npct20:\t\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct20));
    fprintf(fid, 'pct40:\t\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct40));
    fprintf(fid, 'pct1500:\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct1500));
    fprintf(fid, 'pct3000:\t\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.global.Blast.pct3000));
    % Focus
    fprintf(fid, '\n\t\t0-10%%\t10-20%%\t20-100%%\n');
    fprintf(fid, 'Serie Max\t%s\t%s\t%s\n', num2str(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_max),...
                                   num2str(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_max),...
                                   num2str(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_max));
    fprintf(fid, 'Nb de Serie\t%s\t%s\t%s\n', num2str(length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.serie_dur)),...
                                   num2str(length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.serie_dur)),...
                                   num2str(length(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.serie_dur)));
    fprintf(fid, 'Total Time (s)\t%s\t%s\t%s\n', num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur)),...
                                   num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur)),...
                                   num2str(round(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur)));
    fprintf(fid, 'Total Time (%%)\t%s\t%s\t%s\n', num2str(round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m10.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))),...
                                   num2str(round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m20.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))),...
                                   num2str(round(100*(GUI.BLAST_Object.bloc{xi_bloc}.focus.m100.total_dur/GUI.BLAST_Object.bloc{xi_bloc}.duration))));
    
    
    
    fprintf(fid, '\n\n\n');
end




fclose(fid);

function ClicClac(hObj,evnt)

global GUI

[pathTRC, nameTRC, ~] = fileparts(GUI.BLAST_Object.Path_TRC);

fignew = figure('Visible','off'); % Invisible figure

switch hObj.Tag(10:13)
    case 'Glob'
        [filename, pathname] = uiputfile([fullfile(pathTRC, nameTRC) '_Bloc' num2str(GUI.Global.Bouton.choix_bloc_plot.Value) '.pdf'], 'Save Figure as');
        newAxes = copyobj(GUI.Global.fig_Features.Plot_RT.Parent,fignew); % Copy the appropriate axes
    case 'Zoom'
        ob_plt = findobj('Tag', ['seconde_plot_RT_' hObj.Tag(15:end)]);
        [filename, pathname] = uiputfile([fullfile(pathTRC, nameTRC) '_Bloc' num2str(GUI.Global.Bouton.choix_bloc_plot.Value) '_Latency_' num2str(round(ob_plt.Parent.XLim(1))) '_' num2str(round(ob_plt.Parent.XLim(2))) '.pdf'], 'Save Figure as');
        newAxes = copyobj(ob_plt.Parent,fignew); % Copy the appropriate axes
    case 'Bloc'
        [filename, pathname] = uiputfile([fullfile(pathTRC, nameTRC) '_Resum.pdf'], 'Save Figure as');
        GUI.fig_AIC_Pos.main.PaperOrientation = 'landscape';
        fignew = GUI.fig_AIC_Pos.main;
end

%set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
fignew.PaperOrientation = 'landscape';

saveas(fignew, fullfile(pathname, filename))
delete(fignew);

function Rapport(hObj,evnt)


global GUI

[path_rapport, name_rapport, ~] = fileparts(GUI.BLAST_Object.Path_TRC);

len_events = length(GUI.BLAST_Object.bloc{1}.sample.Blast.lat);

[filename, pathname] = uiputfile(fullfile(path_rapport, [name_rapport, '_Rapport_', num2str(len_events) 'events.pdf']), 'Save BLAST Rapport as');

clear len_events


GUI.Global.fig.Children(1).Visible = 'off';
GUI.Global.fig.Children(2).Visible = 'off';
GUI.Global.fig.Children(3).Visible = 'off';
GUI.Global.fig.Children(4).Visible = 'off';
GUI.Global.fig.Children(5).Visible = 'off';
GUI.Global.fig.Children(6).Visible = 'off';
GUI.Global.fig.Children(7).Visible = 'off';
GUI.Global.fig.Children(8).Visible = 'off';
GUI.Global.fig.Children(9).Visible = 'off';
GUI.Global.fig.Children(10).Visible = 'off';
%ax2.Position(3) = ax2.Position(3)*2;


% la on lancera la fonction publish() mais pour l'instant
options = struct('format','pdf',...
                 'outputDir', pathname,...
                 'showCode', false,...
                 'maxWidth', 500);
publish('Buil_Rapport_V02.m', options)

movefile([pathname, 'Buil_Rapport_V02.pdf'], [pathname, filename])

GUI.Global.fig.Children(1).Visible = 'on';
GUI.Global.fig.Children(2).Visible = 'on';
GUI.Global.fig.Children(3).Visible = 'on';
GUI.Global.fig.Children(4).Visible = 'on';
GUI.Global.fig.Children(5).Visible = 'on';
GUI.Global.fig.Children(6).Visible = 'on';
GUI.Global.fig.Children(7).Visible = 'on';
GUI.Global.fig.Children(8).Visible = 'on';
GUI.Global.fig.Children(9).Visible = 'on';
GUI.Global.fig.Children(10).Visible = 'on';

function Launch_RTB(hObj,evnt)
        
global GUI

     
switch GUI.Source
    
    case 'TRC_Blast'
        
        cfg = [];
        cfg.dataset  = GUI.BLAST_Object.Path_TRC;
        cfg.channel  = 'all'; 
        cfg.trl      = [1 2 0];
        data         = ft_preprocessing(cfg);
        GUI.RTB.channel_all  = data.label;
        clear data cfg
    
    case 'MAT_Neurofeedback'
        
        [file_path, file_name, file_ext] = fileparts(GUI.BLAST_Object.Path_TRC);
        load(fullfile(file_path, [file_name(1:end-2), 'OE', file_ext]))

        % There is a difference between control ChannelLabel and TDAH ChannelLabel
        GUI.RTB.channel_all = cellstr(HDR.ChannelLabel);
        clear file_path file_name file_ext HDR DATA
        
    case 'Vamp_Neurofeedback'

        [Fpath,Fname,~] = fileparts([GUI.BLAST_Object.Path_TRC, '_0001']);        

        % only .vhdr are reading by FT
        copyfile(fullfile(Fpath, [Fname,'.ahdr']),...
                 fullfile(Fpath, [Fname,'.vhdr']));
             
        HDR = ft_read_header([GUI.BLAST_Object.Path_TRC, '_0001.vhdr']);
        GUI.RTB.channel_all  = HDR.label;

        clear HDR
        delete(fullfile(Fpath, [Fname,'.vhdr']))
        clear Fpath Fname

end


artefact_param = figure('Position', [50, 100, 500, 500],...
                        'Color', GUI.Colors(1,:), 'MenuBar', 'no',...
                        'Tag', 'artefact_param');

uicontrol('Style', 'text', 'String', {'Select artefact correction'},...
    'Position', [20 465 230 20], 'HorizontalAlignment', 'left',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.methodo = uicontrol('Style', 'popupmenu', 'String', {'Threshold (microVolt)', 'Threshold (%)', 'ICA', 'ICA + Thesh'},...
    'Position', [10 360 200 99], 'FontSize', 15, 'Tag', 'param_band_meth');

% Selection channels to compute RTB
uicontrol('Style', 'text', 'String', {'Channels'},...
    'Position', [300 465 180 20], 'HorizontalAlignment', 'center',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));
uicontrol('Style', 'text', 'String', {'RTB'},...
    'Position', [300 440 75 20], 'HorizontalAlignment', 'right',...
    'FontSize', 12, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.chan_front = uicontrol('Style', 'listbox',...
    'String', intersect(GUI.RTB.channel_all, {'Fz', 'Cz'}), ...
    'Position',[300 365 75 75], 'FontSize', 15);
uicontrol('Style', 'text', 'String', {'All Channels'},...
    'Position', [400 440 75 20], 'HorizontalAlignment', 'right',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));

GUI.RTB.tmp_chx_chan = uicontrol('Style', 'listbox', 'String', GUI.RTB.channel_all,...
    'Position', [400 365 75 75], 'FontSize', 12, "max", 10,...
    'callback', @(~,~)set(GUI.RTB.param_artefact.chan_front,'String', Select_cell(GUI.RTB.channel_all)));


% RTB slide windows features
uicontrol('Style', 'text', 'String', {'RTB Window'},...
    'Position', [300 320 110 20], 'HorizontalAlignment', 'left',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.RTB_win_lengt = uicontrol('Style', 'edit', 'String', '3',...
    'Position', [300 290 50 20], 'FontSize', 15, 'Tag', 'param_RTB_win_size');
uicontrol('Style', 'text', 'String', {'second'},...
    'Position', [355 290 60 20], 'HorizontalAlignment', 'left',...
    'FontSize', 12, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.RTB_win_overlay = uicontrol('Style', 'edit', 'String', '0.875',...
    'Position', [300 260 50 20], 'FontSize', 15, 'Tag', 'param_RTB_win_overlay');
uicontrol('Style', 'text', 'String', {'% Overlay'},...
    'Position', [355 260 60 20], 'HorizontalAlignment', 'left',...
    'FontSize', 12, 'BackgroundColor', GUI.Colors(1,:));


% Artefact detection features
uicontrol('Style', 'frame', 'Position', [1 280 250 130],...
    'BackgroundColor', GUI.Colors(1,:), 'ForegroundColor', [.2 .2 .2]);
uicontrol('Style', 'text', 'String', {'Theshold'},...
    'Position', [10 400 80 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.th_abs_OY = uicontrol('Style', 'edit', 'String', '200',...
    'Position', [20 375 60 20], 'FontSize', 15, 'Tag', 'param_th_abs_OY');
uicontrol('Style', 'text', 'String', {'(micro Volt)'},...
    'Position', [81 375 100 20], 'FontSize', 10,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
uicontrol('Style', 'text', 'String', {'Open Eyes'},...
    'Position', [150 377 100 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.th_abs_CY = uicontrol('Style', 'edit', 'String', '200',...
    'Position', [20 345 60 20], 'FontSize', 15, 'Tag', 'param_th_abs_CY');
uicontrol('Style', 'text', 'String', {'(micro Volt)'},...
    'Position', [81 345 100 20], 'FontSize', 10,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
uicontrol('Style', 'text', 'String', {'Close Eyes'},...
    'Position', [150 347 100 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));

% percent of rejection
uicontrol('Style', 'frame', 'Position', [5 330 220 1],...
    'BackgroundColor', [.5 .5 .5], 'ForegroundColor', [.5 .5 .5]);
GUI.RTB.param_artefact.pct = uicontrol('Style', 'edit', 'String', '0.2',...
    'Position', [20 295 60 20], 'FontSize', 15, 'Tag', 'param_pct');
uicontrol('Style', 'text', 'String', {'(% of signal rejected)'},...
    'Position', [81 295 100 20], 'FontSize', 10,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));

% ICA
GUI.RTB.param_artefact.ICA = {};
uicontrol('Style', 'frame', 'Position', [1 187 250 70],...
    'BackgroundColor', GUI.Colors(1,:), 'ForegroundColor', [.2 .2 .2]);
uicontrol('Style', 'text', 'String', {'ICA Components'},...
    'Position', [10 248 130 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.ICA.Open = uicontrol('Style', 'edit', 'String', 'auto',...
    'Position', [20 220 60 20], 'FontSize', 15, 'Tag', 'ICA Open Eyes');
uicontrol('Style', 'text', 'String', {'Open Eyes'},...
    'Position', [90 222 100 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_artefact.ICA = uicontrol('Style', 'edit', 'String', 'auto',...
    'Position', [20 195 60 20], 'FontSize', 15, 'Tag', 'ICA Close Eyes');
uicontrol('Style', 'text', 'String', {'Close Eyes'},...
    'Position', [90 197 100 20], 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'BackgroundColor', GUI.Colors(1,:));

% Select Spectral correction 
uicontrol('Style', 'frame', 'Position', [1 100 250 50],...
    'BackgroundColor', GUI.Colors(1,:), 'ForegroundColor', [.2 .2 .2]);
uicontrol('Style', 'text', 'String', {'Spectral correction'},...
    'Position', [10 120 140 40], 'HorizontalAlignment', 'left',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.spectral_correction = uicontrol('Style', 'popupmenu', 'String', {'1/f', '1/f + mixt gauss', 'None'},...
    'Position', [10 100 100 30], 'FontSize', 15, 'Tag', 'spectral_correction');



% Select frequency bands definition
uicontrol('Style', 'frame', 'Position', [1 10 250 50],...
    'BackgroundColor', GUI.Colors(1,:), 'ForegroundColor', [.2 .2 .2]);
uicontrol('Style', 'text', 'String', {'Frequency bands definition'},...
    'Position', [10 30 200 40], 'HorizontalAlignment', 'left',...
    'FontSize', 15, 'BackgroundColor', GUI.Colors(1,:));
GUI.RTB.param_band = uicontrol('Style', 'popupmenu', 'String', {'fix', 'Adapt Lansbergen11', 'Adapt Mensia19'},...
    'Position', [10 10 100 30], 'FontSize', 15, 'Tag', 'param_band_meth');

% RTB
uicontrol('Style', 'pushbutton', 'String', 'RTB',...
    'Position', [320 10 60 50], 'FontSize', 15,...
    'BackgroundColor', GUI.Colors(2,:), 'Tag', 'compute_RTB',...
    'callback', 'global GUI, Theta_Beta_Ratio_V05()');
% Save
uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Position', [390 10 60 50], 'FontSize', 15,...
    'BackgroundColor', GUI.Colors(2,:), 'Tag', 'save_RTB',...
    'callback', 'global GUI, Theta_Beta_Ratio_V05()');
% EXIT
uicontrol('Style', 'pushbutton', 'String', 'X',...
    'Position', [460 10 30 30], 'FontSize', 20,...
    'BackgroundColor', [1 0 0], 'Tag', 'quite_RTB',...
    'callback', 'delete(findobj(''Tag'', ''artefact_param''))');

function Extract_Multi_RTB(hObj,evnt)



% Cette fonction va extraire des métrics de la structure de sortie de l'outil de calcule du RTB
% OUTPUT:
%     - nom du sujet
%     
%     Pour chaque enregistrement CY et OY
%     - RTB sur le spectre moyen
%     - RTB sur la moyenne des RTB par fenêtre
%     - Theta sur le spectre moyen
%     - Theta sur la moyenne des RTB par fenêtre
%     - Beta sur le spectre moyen
%     - Beta sur la moyenne des RTB par fenêtre
%     - Alpha sur le spectre moyen
%     - Alpha sur la moyenne des RTB par fenêtre
%     - nombre de fenêtres
%     
%     - peak d'alpha sur le spectre moyen
%     - bande de frequence du theta
%     - bande de frequence du beta 
%     - taille des fennetres
    

% Select all files to import
liste_files = {};
[liste_files.name, liste_files.path, ~] = uigetfile('*_metrics.mat', 'Select output files', 'Multiselect', 'on');

% creat the output file
output_file = {};
[output_file.name, output_file.path, ~] = uiputfile('RTB.txt', 'Save as');

fid = fopen(fullfile(output_file.path, output_file.name), 'w+');
fprintf(fid, 'ID_Suj\tAlpha_Peak\tTheta_Band\tBeta_Band\t');
fprintf(fid, 'OE_Nb_Win\tOE_RTB_Spect_mean\tOE_RTB_Win_Med\tOE_Theta_Spect_mean\tOE_Theta_Win_Med\tOE_Beta_Spect_mean\tOE_Beta_Win_Med\tOE_Alpha_Spect_mean\tOE_Alpha_Win_Med\t');
fprintf(fid, 'CE_Nb_Win\tCE_RTB_Spect_mean\tCE_RTB_Win_Med\tCE_Theta_Spect_mean\tCE_Theta_Win_Med\tCE_Beta_Spect_mean\tCE_Beta_Win_Med\tCE_Alpha_Spect_mean\tCE_Alpha_Win_Med\n');



for xi_file = 1 : numel(liste_files.name)
    
    load(fullfile(liste_files.path, liste_files.name{xi_file}))
    fprintf(fid, '%s\t%s\t%s\t%s\t',...
                liste_files.name{xi_file},...
                num2str(output.alpha_peak),...
                num2str(output.theta_band),...
                num2str(output.beta_band));
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',...
                num2str(length(output.RTB.trial.OY)),...
                num2str(output.RTB.mean.OY),...
                num2str(median(output.RTB.trial.OY)),...
                num2str(output.theta.mean.OY),...
                num2str(median(output.theta.trial.OY)),...
                num2str(output.beta.mean.OY),...
                num2str(median(output.beta.trial.OY)),...
                num2str(output.alpha.mean.OY),...
                num2str(median(output.alpha.trial.OY)));
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
                num2str(length(output.RTB.trial.CY)),...
                num2str(output.RTB.mean.CY),...
                num2str(median(output.RTB.trial.CY)),...
                num2str(output.theta.mean.CY),...
                num2str(median(output.theta.trial.CY)),...
                num2str(output.beta.mean.CY),...
                num2str(median(output.beta.trial.CY)),...
                num2str(output.alpha.mean.CY),...
                num2str(median(output.alpha.trial.CY)));
end, clear xi_file


fclose(fid);



function Change_bloc_display(hObj,evnt)

global GUI

bloc_selected = UI_specify_blocs(length(GUI.BLAST_Object.bloc));

% ajust liste des blocs
ob.cb = findobj('Tag', 'choix_bloc_plot');
ob.cb.String = strsplit(num2str([1:length(bloc_selected)]));

GUI.BLAST_Object.bloc = GUI.BLAST_Object.bloc(str2num(str2mat(bloc_selected)));


Plot_info_refresh()


function output = Select_cell(input_cell)

global GUI 

input_cell = cell2dataset(input_cell, 'ReadVarNames', 0);
output = input_cell.input_cell1(GUI.RTB.tmp_chx_chan.Value);


function cachecache()

%%%%%%%%%%%%
% Connerie pour que le packaging en Appli n'oublie pas cette fonction
%
%%%%%%%%%%%
Buil_Rapport_V02()
Theta_Beta_Ratio_V04()

% ne pas oublier d'ajouter me fichier
% EEG1020.lay
% /Users/romain.bouet/Datas/Sauvegarde/Programmes_Matlab/Fieldtrip/fieldtrip-20200302/template/layout/EEG1020.lay
%
% ne pas oublier de lancer une analyse ICA avant de compiler pour que les
% fichier runica soit dans le path