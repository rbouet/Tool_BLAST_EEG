
function output = UI_specify_items(Bloc_size_min, Bloc_size_max)


global bt_start_first bt_start_last bt_start_specify
global  bt_define_start bt_define_nb_item
global GUI

nb_item_defaut = Bloc_size_max;

GUI.tmp.start = 1;
GUI.tmp.nb_items = nb_item_defaut;
GUI.tmp.size_max = Bloc_size_max;


fig = uifigure('Name', 'BLAST Items Window of interest',...
               'Position',[380 678 398 160],...
               'Color', GUI.Colors(1,:));

bg = uibuttongroup(fig,'Position',[37 53 123 85],...
                   'SelectionChangedFcn',@bselection,...
               'BackgroundColor', GUI.Colors(1,:));
bt_start_first   = uitogglebutton(bg,'Position',[10 50 100 22],...
                                  'text', 'First', 'Tag', 'start_first',...
                                  'BackgroundColor', GUI.Colors(4,:));
bt_start_last    = uitogglebutton(bg,'Position',[10 28 100 22],...
                                  'text', 'Last', 'Tag', 'start_last',...
                                  'BackgroundColor', GUI.Colors(4,:));
bt_start_specify = uitogglebutton(bg,'Position',[10 6 100 22],...
                                  'text', 'Specify', 'Tag', 'start_specify',...
                                  'BackgroundColor', GUI.Colors(4,:));

                              
uilabel(fig,'Text','Nb Items','HorizontalAlignment','left', 'Position', [200 90 50 20]);
bt_define_nb_item = uieditfield(fig, 'Position',[255 90 30 20],...
                               'Value', num2str(GUI.tmp.size_max),...
                               'ValueChangedFcn', @Change_nb_item);

                           uilabel(fig,'Text','Onset','HorizontalAlignment','left', 'Position', [200 53 60 20]);
bt_define_start = uieditfield(fig, 'Position',[255 53 30 20],...
                             'Value', num2str(GUI.tmp.start), 'Enable', 'off',...
                             'ValueChangedFcn', @Change_start);
                         
bt_exit = uiswitch(fig, 'toggle', 'Position',[320 60 30 80],...
                   'Items', {'Select', 'Extract'},...
                   'ValueChangedFcn',@switchMoved);
               
uilabel(fig,'Text',['There are ', num2str(Bloc_size_min), ' items MIN'],...
             'HorizontalAlignment','left', 'Position', [30 25 200 20]);
uilabel(fig,'Text',[num2str(Bloc_size_max), ' items MAX'],...
             'HorizontalAlignment','left', 'Position', [84 10 200 20]);
               
               


uiwait(fig)
output = {};
output.start    = GUI.tmp.start;
output.nb_items = GUI.tmp.nb_items;
GUI = rmfield(GUI, 'tmp');


function bselection(source, event)

global bt_define_start bt_define_nb_item
global GUI

switch source.SelectedObject.Tag
    case 'start_first'
        bt_define_start.Enable = 'off';
        bt_define_start.Value = num2str(1);
    case 'start_last'
        bt_define_start.Enable = 'off';        
        bt_define_start.Value = num2str(GUI.tmp.size_max - str2num(bt_define_nb_item.Value)+1);
    case 'start_specify'
        bt_define_start.Enable = 'on';
end
        
GUI.tmp.start = str2num(bt_define_start.Value);
GUI.tmp.nb_items = str2num(bt_define_nb_item.Value);

function Change_nb_item(source, event)

global bt_define_start bt_define_nb_item
global GUI

GUI.tmp.start = str2num(bt_define_start.Value);
GUI.tmp.nb_items = str2num(bt_define_nb_item.Value);

function Change_start(source, event)

global bt_define_start bt_define_nb_item
global GUI

GUI.tmp.start = str2num(bt_define_start.Value);
GUI.tmp.nb_items = str2num(bt_define_nb_item.Value);

function switchMoved(source, event)

pause(2)
delete(source.Parent)

