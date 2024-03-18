
function output = UI_specify_blocs(nb_bloc_init)

global GUI
global bt_blocs_init bt_blocs_select


fig = uifigure('Name', 'BLAST Blocs',...
               'Position',[380 678 298 200],...
               'Color', GUI.Colors(1,:));

                              
uilabel(fig,'Text','Select blocs to display',...
             'HorizontalAlignment','center',...
             'Position', [1 170 298 20],...
             'FontSize', 15)

uilabel(fig,'Text','Initial',...
             'HorizontalAlignment','left',...
             'Position', [55 125 50 20]);
bt_blocs_init = uilistbox(fig, 'Position',[55 10 50 110],...
                          'Multiselect', 'on',...
                          'Items', strsplit(num2str([1:nb_bloc_init])),...
                          'ValueChangedFcn', @Change_bloc);
                      
uilabel(fig,'Text', 'Selected',...
            'HorizontalAlignment','left',...
            'Position', [110 125 50 20]);
bt_blocs_select = uilistbox(fig, 'Position',[110 10 50 110],...
                          'Multiselect', 'on',...
                          'Items', {''},...
                          'ValueChangedFcn', '');

                      
                      
                         
bt_exit = uiswitch(fig, 'toggle', 'Position',[230 60 30 80],...
                   'Items', {'Select', 'Extract'},...
                   'ValueChangedFcn',@switchMoved);


uiwait(fig)
output = GUI.tmp;
GUI = rmfield(GUI, 'tmp');



function Change_bloc(source, event)

global GUI
global bt_blocs_init bt_blocs_select

bt_blocs_select.Items = bt_blocs_init.Value;


function switchMoved(source, event)

global GUI
global bt_blocs_init bt_blocs_select

GUI.tmp = bt_blocs_select.Items;
pause(2)
delete(source.Parent)

