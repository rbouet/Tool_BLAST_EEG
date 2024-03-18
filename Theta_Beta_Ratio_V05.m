
function Theta_Beta_Ratio()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% meth_artefact : scalar, signal correction method
%                 1 : rejection based on fixed thresold
%                 2 : rejection based on fixed SD
%                 3 : rejection of signal percentage 
%                 4 : ASR (not yet)
%
% meth_RTB : string, band frequency definition
%            'fix'   : fixed frequency band
%            'adapt' : bands were de?ned using alpha as anchor point
%           
% MAJ
%   - 19/09/2021    RB
%               - RTB channels choice
%               - complete RTB, beta, theta and fft save output
%               
%   - 25/01/2022    RB
%               - import neurofeedback Vamp file
%               - create function to create cfg structure for FT structure
%               - remove detrend and add hight pass filter (1Hz)
%
%   - 07/06/2022    RB
%               - Add FOOF spectral correction for RTB 
%               - Add  artefact rejection based on spectrum
%               - Add Mensia's definition of frenquency bands adapted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global GUI

bt_push = gco;


%% Features description

% remove the 5 first secondes
GUI.RTB.param_signal.remove = 5;
% windows lenght, secondes
GUI.RTB.param_signal.win_w = str2num(get(GUI.RTB.param_artefact.RTB_win_lengt, 'String'));
% window overlap
GUI.RTB.param_signal.win_ov = str2num(get(GUI.RTB.param_artefact.RTB_win_overlay, 'String'));
% correlation to detecte ICA components
corr_eog_th    = .4;


% channels of interest
chan_front.label = get(GUI.RTB.param_artefact.chan_front, 'String');
chan_occi.label  = {'P3', 'P4', 'Pz', 'O1', 'O2'};
chan_eog.label   = {'EO+', 'EOGD', 'EOGG'};
%chan_ica.label   = GUI.RTB.channel_all;
%chan_ica.label   = {'Fp1', 'Fp2', 'FC1', 'FC2', 'F7', 'F3', 'F4', 'F8', 'Fz', 'C3', 'Cz', 'C4', 'P3', 'P4', 'Pz', 'O1', 'O2', 'EO+', 'EOGD', 'EOGG'};
chan_ica.label   = {chan_front.label{:}, chan_occi.label{:}, chan_eog.label{:}};


GUI.RTB.output.pct_CY = 1;
GUI.RTB.output.pct_OY = 1;
GUI.RTB.output.ICA_CY = [];
GUI.RTB.output.ICA_OY = [];
        


%% Datas Importation         
switch GUI.Source
    
    case 'TRC_Blast'
        
        % windows of interest
        % oppen eyes
        win_OY = [GUI.BLAST_Object.OY.sample(1),...
            GUI.BLAST_Object.OY.sample(2)];
        % closed eyes
        win_CY = [GUI.BLAST_Object.CY.sample(1),...
            GUI.BLAST_Object.CY.sample(2)];
        
        
        data_ica_OY_continue_raw = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_OY, double(GUI.BLAST_Object.Fs), 0, 'none');
        data_ica_CY_continue_raw = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_CY, double(GUI.BLAST_Object.Fs), 0, 'none');
        data_ica_OY_continue_filt = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_OY, double(GUI.BLAST_Object.Fs), 0, 20);
        data_ica_CY_continue_filt = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_CY, double(GUI.BLAST_Object.Fs), 0, 20);
        
        data_ica_OY_trial = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_OY, double(GUI.BLAST_Object.Fs), 1, 'none');
        data_ica_CY_trial = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_ica, win_CY, double(GUI.BLAST_Object.Fs), 1, 'none');
    
    case 'MAT_Neurofeedback'
        
        data_eog_OY_continue = {};
        data_eog_CY_continue = {};
        
        data_ica_OY_continue_raw = Extract_datas_from_neurofeedback_mat('OE', chan_ica, 0, 'none');
        data_ica_CY_continue_raw = Extract_datas_from_neurofeedback_mat('CE', chan_ica, 0, 'none');
        data_ica_OY_continue_filt = Extract_datas_from_neurofeedback_mat('OE', chan_ica, 0, 20);
        data_ica_CY_continue_filt = Extract_datas_from_neurofeedback_mat('CE', chan_ica, 0, 20);
        
        data_ica_OY_trial = Extract_datas_from_neurofeedback_mat('OE', chan_ica, 1, 'none'); 
        data_ica_CY_trial = Extract_datas_from_neurofeedback_mat('CE', chan_ica, 1, 'none');
   
        
        win_OY = data_ica_OY_continue_raw.sampleinfo+1;
        % closed eyes
        win_CY = data_ica_CY_continue_raw.sampleinfo+1;
        
        
    case 'Vamp_Neurofeedback'
                
        data_eog_OY_continue = {};
        data_eog_CY_continue = {};

        data_ica_OY_continue_raw = Extract_datas_from_neurofeedback_vamp(chan_ica,  'open',  0, 'none');
        data_ica_CY_continue_raw = Extract_datas_from_neurofeedback_vamp(chan_ica,  'close', 0, 'none');
        data_ica_OY_continue_filt = Extract_datas_from_neurofeedback_vamp(chan_ica, 'open',  0, 20);
        data_ica_CY_continue_filt = Extract_datas_from_neurofeedback_vamp(chan_ica, 'close', 0, 20);

        data_ica_OY_trial = Extract_datas_from_neurofeedback_vamp(chan_ica, 'open',  1, 'none');
        data_ica_CY_trial = Extract_datas_from_neurofeedback_vamp(chan_ica, 'close', 1, 'none');
        
        win_OY = [1 data_ica_OY_continue_raw.hdr.nSamples];
        GUI.BLAST_Object.OY.sample = win_OY;
        win_CY = [1 data_ica_CY_continue_raw.hdr.nSamples];
        GUI.BLAST_Object.CY.sample = win_CY;
        

end



%% Preprocesing

% Used later to compute the purcentage of signal rejected later
length_trial_OY = numel(data_ica_OY_trial.trial);
length_trial_CY = numel(data_ica_CY_trial.trial);


% Threshold definition
switch GUI.RTB.param_artefact.methodo.String{GUI.RTB.param_artefact.methodo.Value}
    % define threshold according to methodo selected
    
    case {'Threshold (microVolt)', 'ICA + Thesh'}
        % remove samples upper than threshold
        GUI.RTB.output.threshold_CY = str2num(GUI.RTB.param_artefact.th_abs_CY.String);
        GUI.RTB.output.threshold_OY = str2num(GUI.RTB.param_artefact.th_abs_CY.String);
        fprintf('\nRemove artefact upper than \n\t\t\t\tClosed eyes: %s microVolt\n\t\t\t\tOpened eyes: %s microVolt\n',...
                                        GUI.RTB.param_artefact.th_abs_CY.String,...
                                        GUI.RTB.param_artefact.th_abs_OY.String)
        
        
    case 'Threshold (%)'
        % remove sample to keep a purcentage of signal
        pct_reject_max = str2num(GUI.RTB.param_artefact.pct.String);
        fprintf('\nRemove artefact to keep %s%% of signal\n', num2str(pct_reject_max*100))
        
        flag_CY = 1;
        flag_OY = 1;
        
        for xi_th = 1000 : -5 : 1
            
            to_keep_CY = Reject_seuil_absolu(data_ica_CY_trial, data_ica_CY_continue_raw, xi_th);
            to_keep_OY = Reject_seuil_absolu(data_ica_OY_trial, data_ica_OY_continue_raw, xi_th);
            
            if length(to_keep_CY.trial)/numel(data_ica_CY_trial.trial) < 1-pct_reject_max & flag_CY
                GUI.RTB.output.threshold_CY = xi_th;
                flag_CY = 0;
            end
            if length(to_keep_OY.trial)/numel(data_ica_OY_trial.trial) < 1-pct_reject_max & flag_OY
                GUI.RTB.output.threshold_OY = xi_th;
                flag_OY = 0;
            end
            
        end, clear xi_th flag_OY flag_CY pct_reject_max
        fprintf('\nRemove artefact upper than \n\t\t\t\tClosed eyes: %s microVolt\n\t\t\t\tOpened eyes: %s microVolt\n',...
                                        GUI.RTB.param_artefact.th_abs_CY.String,...
                                        GUI.RTB.param_artefact.th_abs_OY.String)
        
    case 'ICA'
        
        GUI.RTB.output.threshold_CY = 10000000;
        GUI.RTB.output.threshold_OY = 10000000;
end


% find bad trial and bad samples according to threshold
to_keep_CY = Reject_seuil_absolu(data_ica_CY_trial, data_ica_CY_continue_raw, GUI.RTB.output.threshold_CY);
to_keep_OY = Reject_seuil_absolu(data_ica_OY_trial, data_ica_OY_continue_raw, GUI.RTB.output.threshold_OY);
clear th_abs_CY th_abs_OY

% remove bads trials according to amplitude threshold 
cfg = {};
cfg.trials = to_keep_OY.trial;
data_ica_OY_trial = ft_selectdata(cfg, data_ica_OY_trial);
cfg.trials = to_keep_CY.trial;
data_ica_CY_trial = ft_selectdata(cfg, data_ica_CY_trial);
clear cfg


% Remove bad sample along continue datas
% (used only for rejection display)
no_keep = find(ismember(data_ica_OY_continue_raw.time{1}, to_keep_OY.time) == 0);
data_ica_OY_continue_filt.trial{1}(:,no_keep) = nan;

no_keep = find(ismember(data_ica_CY_continue_raw.time{1}, to_keep_CY.time) == 0);
data_ica_CY_continue_filt.trial{1}(:,no_keep) = nan;



    

% ICA correction
if strcmp(GUI.RTB.param_artefact.methodo.String{GUI.RTB.param_artefact.methodo.Value}, 'ICA') || strcmp(GUI.RTB.param_artefact.methodo.String{GUI.RTB.param_artefact.methodo.Value}, 'ICA + Thesh')
        % ICA
        fprintf('ICA correction\n\n')
        % Open eyes
        data_ica_OY_trial = artefact_eog(data_ica_OY_continue_filt,...
                                         data_ica_OY_trial,...
                                         chan_eog,...
                                         corr_eog_th, 'ICA Open Eyes');
        
        
        % closed eyes
         data_ica_CY_trial = artefact_eog(data_ica_CY_continue_filt,...
                                          data_ica_CY_trial,...
                                          chan_eog,...
                                          corr_eog_th, 'ICA Close Eyes');        
                           
end % if meth_artefact

clear data_ica_CY_continue_filt data_ica_CY_continue_raw 
clear data_ica_OY_continue_filt data_ica_OY_continue_raw 

if (isempty(data_ica_CY_trial.trial) || isempty(data_ica_OY_trial.trial))
    errordlg('No datas after thresholding')
    return
end 

% remove EOG channels
cfg         = [];
cfg.channel = setdiff(chan_ica.label, chan_eog.label);
data_ica_OY_trial    = ft_selectdata(cfg, data_ica_OY_trial);
data_ica_CY_trial    = ft_selectdata(cfg, data_ica_CY_trial);


% Select good trial based on spectrim
to_keep_CY = Select_trial_based_spectrum(data_ica_CY_trial);
to_keep_OY = Select_trial_based_spectrum(data_ica_OY_trial);

% remove bads trials according to amplitude threshold 
cfg = {};
cfg.trials = to_keep_OY;
data_ica_OY_trial = ft_selectdata(cfg, data_ica_OY_trial);
cfg.trials = to_keep_CY;
data_ica_CY_trial = ft_selectdata(cfg, data_ica_CY_trial);
clear cfg to_keep_CY to_keep_OY


% compute the purcentage of signal rejected
GUI.RTB.output.pct_CY = numel(data_ica_CY_trial.trial)/length_trial_CY;
GUI.RTB.output.pct_OY = numel(data_ica_OY_trial.trial)/length_trial_OY;

                      
% check correction and artefact removal
chan_disp.label = {chan_front.label{:}, chan_occi.label{:}};
display_preproc_V03(win_CY, win_OY,...
                    chan_disp, chan_eog,...
                    data_ica_CY_trial, data_ica_OY_trial)
clear chan_disp


        
% 1/f correction
[data_ica_CY_mean, data_ica_CY_trial] = Correction_FOOOF(data_ica_CY_trial);
[data_ica_OY_mean, data_ica_OY_trial] = Correction_FOOOF(data_ica_OY_trial);


switch GUI.RTB.param_band.String{GUI.RTB.param_band.Value}
% frequency band definition
        
    case 'fix'
        theta_band = [4 8];
        beta_band  = [13 21];
        alpha_peak = 10;
                
    case 'Adapt Lansbergen11'
        
        [theta_band beta_band alpha_peak] = Search_band_epoked(data_ica_CY_trial, data_ica_OY_trial,...
                                                               chan_occi, 'Lansbergen11');

    case'Adapt Mensia19'
        
        [theta_band beta_band alpha_peak] = Search_band_epoked(data_ica_CY_mean, data_ica_OY_mean,...
                                                               chan_occi, 'Mensia19');

end % switch

GUI.RTB.output.theta_band = theta_band;
GUI.RTB.output.beta_band  = beta_band;
GUI.RTB.output.alpha_peak = alpha_peak;

% frontal labels selection
cfg         = [];
cfg.channel = chan_front.label;
data_ica_OY_trial    = ft_selectdata(cfg, data_ica_OY_trial);
data_ica_CY_trial    = ft_selectdata(cfg, data_ica_CY_trial);
data_ica_OY_mean     = ft_selectdata(cfg, data_ica_OY_mean);
data_ica_CY_mean     = ft_selectdata(cfg, data_ica_CY_mean);

% Comput RTB
[GUI.RTB.output.theta.trial.CY, GUI.RTB.output.beta.trial.CY, GUI.RTB.output.RTB.trial.CY, GUI.RTB.output.alpha.trial.CY] = Compute_RTB(squeeze(nanmean(data_ica_CY_trial.powspctrm, 2)),...
                                                                                                                                      data_ica_CY_trial.freq, theta_band, beta_band, alpha_peak);
[GUI.RTB.output.theta.mean.CY,  GUI.RTB.output.beta.mean.CY,  GUI.RTB.output.RTB.mean.CY,  GUI.RTB.output.alpha.mean.CY]  = Compute_RTB(squeeze(nanmean(data_ica_CY_mean.powspctrm, 1)),...
                                                                                                                                      data_ica_CY_mean.freq,  theta_band, beta_band, alpha_peak);

[GUI.RTB.output.theta.trial.OY, GUI.RTB.output.beta.trial.OY, GUI.RTB.output.RTB.trial.OY, GUI.RTB.output.alpha.trial.OY] = Compute_RTB(squeeze(nanmean(data_ica_OY_trial.powspctrm, 2)),...
                                                                                                                                      data_ica_OY_trial.freq, theta_band, beta_band, alpha_peak);
[GUI.RTB.output.theta.mean.OY,  GUI.RTB.output.beta.mean.OY,  GUI.RTB.output.RTB.mean.OY,  GUI.RTB.output.alpha.mean.OY]  = Compute_RTB(squeeze(nanmean(data_ica_OY_mean.powspctrm, 1)),...
                                                                                                                                      data_ica_OY_mean.freq,  theta_band, beta_band, alpha_peak);

GUI.RTB.output.ff.OY.trial = data_ica_OY_trial;
GUI.RTB.output.ff.CY.trial = data_ica_CY_trial;
GUI.RTB.output.ff.OY.mean = data_ica_OY_mean;
GUI.RTB.output.ff.CY.mean = data_ica_CY_mean;




Display_output()

if strcmp(bt_push.Tag, 'save_RTB')
    Save_output()
    
    return
end

    


function chan = Extract_chan_name_FS(GUI)
% Extraction des labels des electrodes avec les fonctions de FS  
        
EEG_hd = TRC(GUI.BLAST_Object.Path_TRC);


% selection des channels
chan = {};
chan.label = {};
chan.id = [];

for xi_chan = 1 : length(EEG_hd.m_electrodes)
    % si le channel est un capteur EEG
    if EEG_hd.m_electrodes(xi_chan).type == 0
        chan.label{end+1} = EEG_hd.m_electrodes(xi_chan).positivInputLabel;
        chan.id(end+1)    = xi_chan;
    elseif strcmp(EEG_hd.m_electrodes(xi_chan).positivInputLabel, 'EO+')
        chan.label{end+1} = EEG_hd.m_electrodes(xi_chan).positivInputLabel;
        chan.id(end+1)    = xi_chan;
    end% if 
end, clear xi_chan

function data = Extract_datas_Clinic_TRC(filename, chan, win, Fs, epoked, lowpass)
% Epoked datas extraction

global GUI

if epoked
    trl = Define_TRL_epoked(win);
else
    trl = [win(1)+(GUI.RTB.param_signal.remove*Fs),...
           win(2),...
           0];
end % if epoked
   

cfg = [];
cfg.dataset  = filename;
cfg.channel  = chan.label; %'eeg';
% cfg.hpfilter = 'yes';
% cfg.hpfreq   = 1;
if isnumeric(lowpass)
   cfg.lpfilter = 'yes';
   cfg.lpfreq = lowpass;
else
    cfg.lpfilter = 'no';
end
cfg.trl      = trl;
data         = ft_preprocessing(cfg);
data.sampleinfo = data.sampleinfo - data.sampleinfo(1,1);

function data = Extract_datas_from_neurofeedback_vamp(label, eyes, epoked, lowpass)

global GUI

switch eyes
    case 'open'
        [Fpath,Fname,~] = fileparts([GUI.BLAST_Object.Path_TRC, '_0001']);
        
    case 'close'
        [Fpath,Fname,~] = fileparts([GUI.BLAST_Object.Path_TRC, '_0002']);
end


% only .vmrk and .vhdr are reading by FT
copyfile(fullfile(Fpath, [Fname,'.amrk']),...
    fullfile(Fpath, [Fname,'.vmrk']));
copyfile(fullfile(Fpath, [Fname,'.ahdr']),...
    fullfile(Fpath, [Fname,'.vhdr']));

data = {};
% Build Header
data.hdr = ft_read_header(fullfile(Fpath, [Fname,'.vhdr']));
HDR_tmp = data.hdr;
HDR_tmp.nChans = data.hdr.nChans + 1;
HDR_tmp.orig.NumberOfChannels = HDR_tmp.nChans;
HDR_tmp.nSamples = data.hdr.nSamples*data.hdr.nChans/(data.hdr.nChans + 1);
HDR_tmp.orig.nSamples = HDR_tmp.nSamples;


data.label      = data.hdr.label;
data.fsample    = data.hdr.Fs;

% Import datas
data.trial{1} = ft_read_data(fullfile(Fpath, [Fname,'.eeg']),...
                             'header',    HDR_tmp,...
                             'begsample', 1,...
                             'endsample', HDR_tmp.nSamples,...
                             'chanindx',  1:data.hdr.nChans);
clear HDR_tmp

data.sampleinfo = [0 size(data.trial{1}, 2)-1];
data.hdr.nSamples = size(data.trial{1}, 2);
data.time{1}    = 0:1/data.hdr.Fs:(data.hdr.nSamples-1)/data.hdr.Fs;

% HDR.nSamples = size(DATA,1);
% HDR.orig.nSamples = size(DATA,1);

% create cfg, filter and epoke
data = Creat_cfg_for_new_ft_structure(data, label, epoked, lowpass);

delete(fullfile(Fpath, [Fname, '.vhdr']));
delete(fullfile(Fpath, [Fname, '.vmrk']));

function data = Extract_datas_from_neurofeedback_mat(eyes, label, epoked, lowpass)

% file_path doit contenir le nom du fichier SANS le suffixe CE ou OE
% eyes : 'CE' ou 'OE'
% label : label to export     label = chan_front
% epoked : 0 ou 1 si on veut des données epoké ou non
% lowpass : frequence du filtre lowpass: 'none' ou scalar

global GUI

% load MAT file according to eyes
[file_path, file_name, file_ext] = fileparts(GUI.BLAST_Object.Path_TRC);
load(fullfile(file_path, [file_name(1:end-2), eyes, file_ext]))

% There is a difference between control ChannelLabel and TDAH ChannelLabel
HDR.ChannelLabel = cellstr(HDR.ChannelLabel);

% center datas
DATA = DATA-repmat(mean(DATA), size(DATA, 1), 1);


GUI.BLAST_Object.Suj_name = file_name(1:end-3);


data = {};

data.hdr = {};
data.hdr.Fs          = HDR.SamplingFrequency;
data.hdr.nChans      = numel(HDR.ChannelLabel);
data.hdr.nSamples    = HDR.NumberOfSamples;
data.hdr.nSamplesPre = 0;
data.hdr.nTrials     = 1;
data.hdr.label       = HDR.ChannelLabel;
data.hdr.chanunit    = {};
data.hdr.subjectname = GUI.BLAST_Object.Suj_name;
data.hdr.chantype    = {};

for xi = 1 : data.hdr.nChans
    % transformation de Volt a microvolt
    if mean(DATA(:))+std(DATA(:)) < 1
        DATA = DATA.*1e6;
    end
    data.hdr.chanunit{end+1} = 'uV';
    data.hdr.chantype{end+1} = 'eeg';
end, clear xi


data.label = HDR.ChannelLabel;
data.time = {};
data.trial = {};
data.fsample = HDR.SamplingFrequency;
data.sampleinfo = {};

%data.time{1} = [0 : 1/data.fsample : (data.hdr.nSamples-1)/data.fsample];
data.time{1} = double([0 : data.hdr.nSamples-1])/data.fsample;
data.sampleinfo = [0 data.hdr.nSamples-1];
% 
% % find channels within HDR.ChannelLabel
% index = cellfun(@(k) find(strcmp(k, HDR.ChannelLabel)), label.label, 'UniformOutput', false);
% index = cell2mat(index');
% data.trial{1} = DATA(:,index)';
 data.trial{1} = DATA';

data = Creat_cfg_for_new_ft_structure(data, label, epoked, lowpass);

function data = Creat_cfg_for_new_ft_structure(data, label, epoked, lowpass)
        
global GUI

data.cfg = {};
data.cfg.dataset = GUI.BLAST_Object.Path_TRC;
data.cfg.channel = data.label;
data.cfg.lpfilter = 'no';
data.cfg.lpfreq = [];
data.cfg.trl = [1 data.hdr.nSamples 0];
data.cfg.checkpath = 'pedantic';
data.cfg.outputfilepresent=  'overwrite';
data.cfg.toolbox.images = 'compat';
data.cfg.toolbox.stats = 'compat';
data.cfg.toolbox.signal = 'compat';
data.cfg.toolbox.cleanup = {};
data.cfg.progress.noerase = 0;
data.cfg.callinfo.usercfg.dataset = '/Volumes/Backup DD/From_BAIE/Vania/Analyse_Clinic/Stabilo/datas_raw/Cognit-AIC_027DR190318.TRC';
data.cfg.callinfo.usercfg.channel = data.label;
data.cfg.callinfo.usercfg.lpfilter = 'no';
data.cfg.callinfo.usercfg.lpfreq = [];
data.cfg.callinfo.usercfg.trl = [1 data.hdr.nSamples 0];
data.cfg.callinfo.usercfg.checkpath = 'pedantic';
data.cfg.callinfo.usercfg.outputfilepresent = 'overwrite';
data.cfg.callinfo.usercfg.toolbox.images= 'compat';
data.cfg.callinfo.usercfg.toolbox.stats = 'compat';
data.cfg.callinfo.usercfg.toolbox.signal = 'compat';
data.cfg.callinfo.usercfg.toolbox.cleanup = {};
data.cfg.callinfo.usercfg.progress.noerase = 0;
data.cfg.callinfo.fieldtrip = '20200302';
data.cfg.callinfo.matlab = '9.6.0.1072779 (R2019a)';
data.cfg.callinfo.computer = 'maci64';
data.cfg.callinfo.hostname = 'mac-pro-de-romainbouet';
data.cfg.callinfo.user = 'romain.bouet';
data.cfg.callinfo.pwd = '/Volumes/Backup DD/From_BAIE/Vania/Analyse_Clinic/Stabilo/Scripts';
data.cfg.callinfo.calltime = [2021 3 17 11 26 51.4132];
data.cfg.callinfo.proctime = 0.7909;
data.cfg.callinfo.procmem = 32976896;
data.cfg.version.name = '/Users/romain.bouet/Datas/Sauvegarde/Programmes_Matlab/Fieldtrip/fieldtrip-20200302/ft_preprocessing.m';
data.cfg.version.id = '$Id$';
data.cfg.method = 'trial';
data.cfg.removemcg = 'no';
data.cfg.removeeog = 'no';
data.cfg.precision = 'double';
data.cfg.padding = 0;
data.cfg.paddir = 'both';
data.cfg.headerformat = 'brainvision_vhdr';
data.cfg.dataformat = 'brainvision_eeg';
data.cfg.coordsys = 'head';
data.cfg.coilaccuracy = [];
data.cfg.checkmaxfilter = [];
data.cfg.montage = 'no';
data.cfg.updatesens = 'no';
data.cfg.chantype = {};
data.cfg.dftfilter = 'no';
data.cfg.hpfilter = 'no';
data.cfg.bpfilter = 'no';
data.cfg.bsfilter = 'no';
data.cfg.medianfilter = 'no';
data.cfg.padtype = 'data';
data.cfg.reref = 'no';
data.cfg.refchannel = {};
data.cfg.refmethod = 'avg';
data.cfg.implicitref = [];
data.cfg.feedback = 'text';
data.cfg.datafile = GUI.BLAST_Object.Path_TRC;
data.cfg.headerfile = GUI.BLAST_Object.Path_TRC;
data.cfg.continuous = 'yes';
data.cfg.polyremoval = 'no';
data.cfg.polyorder = 2;
data.cfg.detrend = 'no';
data.cfg.demean = 'no';
data.cfg.baselinewindow = 'all';
data.cfg.lpfiltord = [];
data.cfg.hpfiltord = [];
data.cfg.bpfiltord = [];
data.cfg.bsfiltord = [];
data.cfg.lpfilttype = 'but';
data.cfg.hpfilttype = 'but';
data.cfg.bpfilttype = 'but';
data.cfg.bsfilttype = 'but';
data.cfg.lpfiltdir = 'twopass';
data.cfg.hpfiltdir = 'twopass';
data.cfg.bpfiltdir = 'twopass';
data.cfg.bsfiltdir = 'twopass';
data.cfg.lpinstabilityfix = 'no';
data.cfg.hpinstabilityfix = 'no';
data.cfg.bpinstabilityfix = 'no';
data.cfg.bsinstabilityfix = 'no';
data.cfg.lpfiltdf = [];
data.cfg.hpfiltdf = [];
data.cfg.bpfiltdf = [];
data.cfg.bsfiltdf = [];
data.cfg.lpfiltwintype = 'hamming';
data.cfg.hpfiltwintype = 'hamming';
data.cfg.bpfiltwintype = 'hamming';
data.cfg.bsfiltwintype = 'hamming';
data.cfg.lpfiltdev = [];
data.cfg.hpfiltdev = [];
data.cfg.bpfiltdev = [];
data.cfg.bsfiltdev = [];
data.cfg.plotfiltresp = 'no';
data.cfg.usefftfilt = 'no';
data.cfg.medianfiltord = 9;
data.cfg.dftfreq = [50 100 150];
data.cfg.hilbert = 'no';
data.cfg.derivative = 'no';
data.cfg.rectify = 'no';
data.cfg.boxcar = 'no';
data.cfg.absdiff = 'no';
data.cfg.conv = 'no';
data.cfg.dftinvert = 'no';
data.cfg.standardize = 'no';
data.cfg.denoise = '';
data.cfg.subspace = [];
data.cfg.custom = '';
data.cfg.resample = '';
data.cfg.previous = {};

% channel selection
cfg = {};
cfg.channel = label.label;
data = ft_selectdata(cfg, data);
clear cfg

% cfg = [];
% cfg.detrend = 'yes';
% data = ft_preprocessing(cfg, data);
% clear cfg    

cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data = ft_preprocessing(cfg, data);
clear cfg

if isnumeric(lowpass)
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = lowpass;
    data = ft_preprocessing(cfg, data);
    clear cfg
end

if epoked == 1
    cfg = {};
    cfg.trl = Define_TRL_epoked(data.sampleinfo);

    data = ft_redefinetrial(cfg, data);
    clear cfg
end

function trl = Define_TRL_epoked(win)
% Construction of the EPOKED trl matrix

global GUI

Fs     = double(GUI.BLAST_Object.Fs);
remove = GUI.RTB.param_signal.remove;
win_w  = GUI.RTB.param_signal.win_w;
win_ov = GUI.RTB.param_signal.win_ov;

start_win = [win(1)+round(remove*Fs) : round((win_w*Fs)*(1-win_ov)) : win(2)-round(win_w*Fs)];
stop_win = start_win + round(win_w*Fs);
trl = [start_win',...
       stop_win',...
       zeros(length(start_win),1)];
    
function data_e = Epoked_continue_data(data)

global GUI

win = [data.time{1}(1) data.time{1}(end)];

win_w  = GUI.RTB.param_signal.win_w;
win_ov = GUI.RTB.param_signal.win_ov;
    
start_win = [win(1) : round(win_w*(1-win_ov)) : round(win(2)-win_w)];
stop_win = start_win + win_w;
trl = [start_win',...
       stop_win',...
       zeros(length(start_win),1)];
    
   
data_e = data;
data_e.trial = {};
data_e.time = {};
data_e.sampleinfo = [];

for xi_trl = 1 : size(trl, 1)
    data_e.trial{xi_trl}           = data.trial{1}(find(data.time{1} == trl(xi_trl,1)):find(data.time{1} == trl(xi_trl,2)));
    data_e.time{xi_trl}            = [0 : 1/data.hdr.Fs : win_w];
    data_e.sampleinfo(xi_trl, 1:2) = [trl(xi_trl, 1:2)]*data.hdr.Fs;
end, clear xi_trl
 
function to_keep = Reject_seuil_absolu(data_trial, data_continue, th)
% if only one sample is upper than TH 
% then the trial is rejected

global GUI

ref_time = 0;
to_keep = {};
to_keep.trial = [];
to_keep.time = [];
to_keep.starts = [];

for xi_trial = 1 : numel(data_trial.trial)
    
      % Signal centering
      data_trial.trial{xi_trial} = data_trial.trial{xi_trial} - mean(data_trial.trial{xi_trial},2);
      
    if max(abs(data_trial.trial{xi_trial}(:))) < th
        to_keep.trial = [to_keep.trial, xi_trial];
        to_keep.time = [to_keep.time data_trial.time{xi_trial}+ref_time];
        to_keep.starts = [to_keep.starts ref_time];
    end
    
    ref_time = ref_time + data_trial.time{1}(end)*(1-GUI.RTB.param_signal.win_ov);

end, clear xi_trial

function th_abs = Sd_to_Th(data, st_th)

% th_abs : median + (st_th * std)
% if there are several channel, we keep only the smaller th_abs


th_abs = min(nanmedian(data.trial{1}, 2) + (nanstd(data.trial{1}, 0, 2)*st_th));

function data_clean = Remove_trial(data, rm_id)

data_clean = data;
data_clean.trial = {};
data_clean.time  = {};
data_clean.sampleinfo = [];

for xi_trial = rm_id
    
    data_clean.trial{end+1} = data.trial{xi_trial};
    data_clean.time{end+1}  = data.time{xi_trial};
    data_clean.sampleinfo(end+1,1:2) = data.sampleinfo(xi_trial, 1:2);
    
end, clear xi_trial

function data_reduc = Dimension_reduc(data)

data_reduc = data;
data_reduc.trial = {};
data_reduc.label = {};

for xi_trial = 1 : numel(data.trial)
    
    data_reduc.trial{xi_trial} = nanmean(data.trial{xi_trial},1);
    data_reduc.label{xi_trial} = strcat(data.label{:});
    
end, clear xi_trial

RTB = [];

function [ff, theta, beta, RTB] = Ratio_by_trial(data, theta_band, beta_band)

% This function 
% - compute FFT
% - extract RTB
% By trial
%
% data : ft structure (signal)
%
% OBSOLETE

ff    = nan(data.fsample, numel(data.trial));
theta = nan(1, numel(data.trial));
beta  = nan(1, numel(data.trial));
RTB   = nan(1, numel(data.trial));

for xi_trial = 1 : numel(data.trial)
    
    ff(:, xi_trial) = abs(fft(data.trial{xi_trial}, data.fsample));
    theta(xi_trial) = sum(log(ff(theta_band(1):theta_band(2), xi_trial)));
    beta(xi_trial)  = sum(log(ff(beta_band(1):beta_band(1), xi_trial)));
    RTB(xi_trial)   = theta(xi_trial)/beta(xi_trial);
        
end, clear xi_trial
        
function [ff, theta, beta, RTB] = Ratio_by_trial_corrected_1f(data, theta_band, beta_band)
% This function compute de RTB based on fft by trial
% Here we don't use the log transformation 
% We fit an spectral 1/f curve to remove the frequency unbalance
%
% OBSOLETE

ff    = nan(data.fsample, numel(data.trial));
theta = nan(1, numel(data.trial));
beta  = nan(1, numel(data.trial));
RTB   = nan(1, numel(data.trial));

% initialisation for correction 1/f
freq_to_fit = [2:4 100:250];    
freq_all    = 1:size(ff,1);
    
for xi_trial = 1 : numel(data.trial)
    
    % spectre
    ff(:, xi_trial) = abs(fft(data.trial{xi_trial}, data.fsample));
    
    % Correction
    % WARNING: correction is remove above 512 hz
    warning off
    fout = fit(freq_to_fit', ff(freq_to_fit, xi_trial), 'a + b/x');
    warning on
    ff_fit = fout.a + fout.b./freq_all;
    ff(:, xi_trial) = ff(:, xi_trial) - ff_fit';
    clear fout ff_fit 
    
    % RTB
    theta(xi_trial) = sum(ff(theta_band(1):theta_band(2), xi_trial));
    beta(xi_trial)  = sum(ff(beta_band(1):beta_band(1), xi_trial));
    RTB(xi_trial)   = theta(xi_trial)/beta(xi_trial);
        
end, clear xi_trial

function [good_one] = Select_trial_based_spectrum(data)

% This function slect the good trials 
% We remove Each trial where spectrum of 1 channel is up to mean+2SD 
%
% Input:
%       data       : FT structure
%       theta_band : vector, [min max]Hz
%       beta_band  : vector, [min max]Hz
%
% Output:
%       good_one : vecto, ID of good trials

% Compute mean SPECTRUM raw data
cfg               = [];
cfg.foilim        = [1 40];
cfg.tapsmofrq     = 1;
cfg.method        = 'mtmfft';
cfg.output        = 'pow';
spc_mean = ft_freqanalysis(cfg, data);

% Compute SPECTRUM trial by trial
cfg.keeptrials    = 'yes'; 
spc_trial = ft_freqanalysis(cfg, data);
clear cfg

% figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 4)

% find trials out of spectral 2sd
m_trial  = squeeze(nanmean(spc_trial.powspctrm, 1));
sd_trial = squeeze(nanstd(spc_trial.powspctrm));
% Define how many outliers samples hz need to reject
nb_sample_hz = (nearest(spc_trial.freq, 2) - nearest(spc_trial.freq, 1)) * 3;    % reject 3HZ outliers

% good_one = find(sum((squeeze(spc_trial.powspctrm) - repmat(m_trial+ (2*sd_trial), size(squeeze(spc_trial.powspctrm),1), 1)) > 0, 2) < nb_sample_hz);

good_one = find(max(sum((spc_trial.powspctrm - permute(repmat(m_trial+ (2*sd_trial),...
                                                              1, 1, size(spc_trial.powspctrm,...
                                                                         1)),...
                                                       [3,1,2])) > 0, 3),...
                    [], 2) < nb_sample_hz);

clear nb_sample_hz
% figure, imagesc((squeeze(spc_trial.powspctrm) - repmat(m_trial+ (2*sd_trial), size(squeeze(spc_trial.powspctrm),1), 1)) > 0)
% figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm(good_one,:,:)))
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 4)

function [spc_mean_correct, spc_trial_correct] = Correction_FOOOF(data)
% This function use FOOOF method to remove 1/f trend and gaussian mixture
% soft : https://github.com/fooof-tools
% ref : https://www.nature.com/articles/s41593-020-00744-x
%
% Here we use the implemented fieldtrip method 
% https://www.fieldtriptoolbox.org/example/fooof/
%
% Input:
%       data       : FT structure
%
% Output:
%       spc_trial_correct      : FT structure, powspectrum trial by trial, 
%                                        FOOF corrected
%       good_one               : vector, trial to keep, id
%       theta_m, beta_m, RTB_m : scalar, theta, beta and RTB 
%                                        compute on the mean
%       theta_t, beta_t, RTB_t : scalar, theta, beta and RTB 
%                                        compute trial by trial


global GUI

% Compute mean SPECTRUM on selected data
cfg               = [];
cfg.foilim        = [1 40];
%cfg.pad           = 'maxperlen';
cfg.tapsmofrq     = 2;
%cfg.taper         = 'dpss';
cfg.method        = 'mtmfft';
cfg.output        = 'pow';
% figure, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 2)
spc_mean = ft_freqanalysis(cfg, data);
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'r', 'LineWidth', 2)

% Compute SPECTRUM trial by trial
cfg.keeptrials    = 'yes'; 
spc_trial = ft_freqanalysis(cfg, data);

%%%%%%%%%% FOOOF %%%%%%%%%%%%%%%%%%%%%
fprintf('\nFOOF spectral correction:\t')
switch GUI.RTB.spectral_correction.String{GUI.RTB.spectral_correction.Value}
    case '1/f'
        % compute 1/f trend
        fprintf('1/f\n\n')
        cfg.output        = 'fooof_aperiodic';     
        cfg.keeptrials    = 'no'; 
        spc_component = ft_freqanalysis(cfg, data);
        % figure, plot(trend_1f.freq, trend_1f.powspctrm, 'k', 'LineWidth', 2)
        
        % Correct the mean 
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        spc_mean_correct  = ft_math(cfg, spc_component, spc_mean);
        clear cfg

        % Correct the trial by trial spectrum
        spc_trial_correct = spc_trial;
        spc_trial_correct.powspctrm = spc_trial.powspctrm - permute(repmat(spc_component.powspctrm, 1, 1, size(spc_trial.powspctrm,1)),...
                                                                    [3,1,2]);
        
        
        % figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
        % figure, plot(spc_trial_correct.freq, spc_trial_correct.powspctrm)
        clear spc_component

        
    case '1/f + mixt gauss'
        % estimation of 1/f trend and mixture of gaussians
        fprintf('1/f + Gaussian Mixture\n\n')
        cfg.output        = 'fooof';
        cfg.keeptrials    = 'no'; 
        spc_component = ft_freqanalysis(cfg, data);     
        clear cfg
        % figure, plot(spc_all_component.freq, spc_all_component.powspctrm, 'k', 'LineWidth', 2)
        
        % Correct the mean 
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        spc_mean_correct  = ft_math(cfg, spc_component, spc_mean);
        clear cfg

        % Correct the trial by trial spectrum
        spc_trial_correct = spc_trial;
        spc_trial_correct.powspctrm = spc_trial.powspctrm - permute(repmat(spc_component.powspctrm, 1, 1, size(spc_trial.powspctrm,1)),...
                                                                    [3,1,2]);
        % figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
        % figure, plot(spc_trial_correct.freq, spc_trial_correct.powspctrm)
        clear spc_component
        
    case 'None'
        % estimation of 1/f trend and mixture of gaussians
        fprintf('None\n\n')
        spc_mean_correct  = spc_mean;
        spc_trial_correct = spc_trial;
        
end

function [spc_trial_correct, good_one,...
          theta_m, beta_m, RTB_m,...
          theta_t, beta_t, RTB_t,...
          aplha_m, aplha_t] = Ratio_by_trial_corrected_FOOOF(data, theta_band, beta_band, alpha_peak)
% This function use FOOOF method to remove 1/f trend and gaussian mixture
% soft : https://github.com/fooof-tools
% ref : https://www.nature.com/articles/s41593-020-00744-x
%
% Here we use the implemented fieldtrip method 
% https://www.fieldtriptoolbox.org/example/fooof/
%
% Input:
%       data       : FT structure
%       theta_band : vector, [min max]Hz
%       beta_band  : vector, [min max]Hz
%
% Output:
%       spc_trial_correct      : FT structure, powspectrum trial by trial, 
%                                        FOOF corrected
%       good_one               : vector, trial to keep, id
%       theta_m, beta_m, RTB_m : scalar, theta, beta and RTB 
%                                        compute on the mean
%       theta_t, beta_t, RTB_t : scalar, theta, beta and RTB 
%                                        compute trial by trial


global GUI

% There is an error with new fieldtrip function
% we have to modify data structure
% truc = data.label{1};
% data.label = {};
% data.label{1} = truc;
% clear truc


% Compute mean SPECTRUM raw data
cfg               = [];
cfg.foilim        = [1 40];
cfg.tapsmofrq     = 3;
cfg.method        = 'mtmfft';
cfg.output        = 'pow';
spc_mean = ft_freqanalysis(cfg, data);

% Compute SPECTRUM trial by trial
cfg.keeptrials    = 'yes'; 
spc_trial = ft_freqanalysis(cfg, data);
clear cfg

% figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 4)

% find trials out of spectral 2sd
m_trial  = squeeze(nanmean(spc_trial.powspctrm, 1));
sd_trial = squeeze(nanstd(spc_trial.powspctrm));
% Define how many outliers samples hz need to reject
nb_sample_hz = (nearest(spc_trial.freq, 2) - nearest(spc_trial.freq, 1)) * 3;    % reject 3HZ outliers

% good_one = find(sum((squeeze(spc_trial.powspctrm) - repmat(m_trial+ (2*sd_trial), size(squeeze(spc_trial.powspctrm),1), 1)) > 0, 2) < nb_sample_hz);

good_one = find(max(sum((spc_trial.powspctrm - permute(repmat(m_trial+ (2*sd_trial),...
                                                              1, 1, size(spc_trial.powspctrm,...
                                                                         1)),...
                                                       [3,1,2])) > 0, 3),...
                    [], 2) < nb_sample_hz);

clear nb_sample_hz
% figure, imagesc((squeeze(spc_trial.powspctrm) - repmat(m_trial+ (2*sd_trial), size(squeeze(spc_trial.powspctrm),1), 1)) > 0)
% figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm(good_one,:,:)))
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 4)



% Compute mean SPECTRUM on selected data
cfg               = [];
cfg.foilim        = [1 40];
cfg.tapsmofrq     = 1;
cfg.method        = 'mtmfft';
cfg.output        = 'pow';
cfg.trials         = good_one;
% figure, plot(spc_mean.freq, spc_mean.powspctrm, 'k', 'LineWidth', 2)
spc_mean = ft_freqanalysis(cfg, data);
% hold on, plot(spc_mean.freq, spc_mean.powspctrm, 'r', 'LineWidth', 2)

% Compute SPECTRUM trial by trial
cfg.keeptrials    = 'yes'; 
cfg.trials         = good_one;
spc_trial = ft_freqanalysis(cfg, data);

%%%%%%%%%% FOOOF %%%%%%%%%%%%%%%%%%%%%
fprintf('\nFOOF spectral correction:\t')
switch GUI.RTB.spectral_correction.String{GUI.RTB.spectral_correction.Value}
    case '1/f'
        % compute 1/f trend
        fprintf('1/f\n\n')
        cfg.output        = 'fooof_aperiodic';     
        cfg.keeptrials    = 'no'; 
        spc_component = ft_freqanalysis(cfg, data);
        % figure, plot(trend_1f.freq, trend_1f.powspctrm, 'k', 'LineWidth', 2)
        
        % Correct the mean 
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        spc_mean_correct  = ft_math(cfg, spc_component, spc_mean);
        clear cfg

        % Correct the trial by trial spectrum
        spc_trial_correct = spc_trial;
        spc_trial_correct.powspctrm = squeeze(spc_trial.powspctrm) - repmat(spc_component.powspctrm, size(spc_trial.powspctrm,1), 1);
        % figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
        % figure, plot(spc_trial_correct.freq, spc_trial_correct.powspctrm)
        clear spc_component

        
    case '1/f + mixt gauss'
        % estimation of 1/f trend and mixture of gaussians
        fprintf('1/f + Gaussian Mixture\n\n')
        cfg.output        = 'fooof';
        cfg.keeptrials    = 'no'; 
        spc_component = ft_freqanalysis(cfg, data);     
        clear cfg
        % figure, plot(spc_all_component.freq, spc_all_component.powspctrm, 'k', 'LineWidth', 2)
        
        % Correct the mean 
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        spc_mean_correct  = ft_math(cfg, spc_component, spc_mean);
        clear cfg

        % Correct the trial by trial spectrum
        spc_trial_correct = spc_trial;
        spc_trial_correct.powspctrm = squeeze(spc_trial.powspctrm) - repmat(spc_component.powspctrm, size(spc_trial.powspctrm,1), 1);
        % figure, plot(spc_trial.freq, squeeze(spc_trial.powspctrm))
        % figure, plot(spc_trial_correct.freq, spc_trial_correct.powspctrm)
        clear spc_component
        
    case 'None'
        % estimation of 1/f trend and mixture of gaussians
        fprintf('None\n\n')
        spc_mean_correct  = spc_mean;
        spc_trial_correct = spc_trial;
        
end



% % Extract RTB mean
[theta_m, beta_m, RTB_m, aplha_m] = Compute_RTB(spc_mean_correct.powspctrm, spc_mean_correct.freq, theta_band, beta_band, alpha_peak);
[theta_t, beta_t, RTB_t, aplha_t] = Compute_RTB(spc_trial_correct.powspctrm, spc_trial_correct.freq, theta_band, beta_band, alpha_peak);




function [theta, beta, RTB, aplha] = Compute_RTB(spc, freq, theta_band, beta_band, alpha_peak)
% This function extract RTB
% 
% spc   : vector or matrix (trial X freq), spectral power
% freq  : vector, spectral frequency  
% theta : vector, [min max]
% beta  : vector, [min max]

theta = median(spc(:,...
                 nearest(freq, theta_band(1)):nearest(freq, theta_band(2))),...
             2);
beta  = median(spc(:,...
                 nearest(freq, beta_band(1)):nearest(freq, beta_band(2))),...
             2);
RTB   = theta./beta;
aplha = median(spc(:,...
                 nearest(freq, alpha_peak-1):nearest(freq, alpha_peak+1)),...
             2);


function [theta_band beta_band alpha_peak] = Search_band(data_occi_CY, data_occi_OY)
 
% Search alpha peak between 5 and 15 Hz
% bands were de?ned using IAF as anchor point
% theta and beta according to Lansbergen et al 2011


global GUI

data_occi_CY = nanmean(data_occi_CY.trial{1},1);
data_occi_OY = nanmean(data_occi_OY.trial{1},1);

fft_CY = abs(fft(data_occi_CY, GUI.BLAST_Object.Fs));
fft_OY = abs(fft(data_occi_OY, GUI.BLAST_Object.Fs));
fft_diff = fft_CY - fft_OY;

alpha_peak = 4 + find(fft_diff(5:15) == max(fft_diff(5:15)));
theta_band = [0.4*alpha_peak 0.6*alpha_peak];
beta_band = [1.2*alpha_peak 25];

x = [1:25];
figure('Position', [100,100, 500 500],  'MenuBar', 'no', 'Color', GUI.Colors(1,:),...
       'Tag', 'Frequency_power');
plot(x,fft_CY(x), '--b')
hold on, plot(x,fft_OY(x), '--k')
plot(x,fft_diff(x), 'r')
legend({'Eyes Close', 'Eyes Open', 'Difference'})
xlabel('Hz')
      
function [theta_band beta_band alpha_peak] = Search_band_epoked(data_occi_CY, data_occi_OY, chan_occi, method)
 
% Search alpha peak between 5 and 15 Hz
% bands were de?ned using IAF as anchor point
% theta and beta according to Lansbergen et al 2011


global GUI

% EOG channels selection
cfg          = [];
cfg.channel  = chan_occi.label;
data_occi_CY = ft_selectdata(cfg, data_occi_CY);
data_occi_OY = ft_selectdata(cfg, data_occi_OY);
clear cfg


fft_CY = mean(data_occi_CY.powspctrm,1);
fft_OY = mean(data_occi_OY.powspctrm,1);

fft_diff = fft_CY - fft_OY;
freq_id = nearest(data_occi_CY.freq, 7): nearest(data_occi_CY.freq, 12);

alpha_peak = data_occi_CY.freq(freq_id(1)-1 + find(fft_diff(freq_id) == max(fft_diff(freq_id))));

switch method
    case 'Lansbergen11'
        % according to Lansbergen et al 2011
        theta_band = [0.4*alpha_peak 0.8*alpha_peak];
        beta_band = [1.2*alpha_peak 25];
    case 'Mensia19'
        % according to Mensia 2019
        theta_band = [alpha_peak-5 alpha_peak-1];
        beta_band = [alpha_peak+3 alpha_peak+12];
        
end



function output = slide_ratio(data_eeg, artifact)

% Compute spctrum on 2s signal with 1s slide windows

slide_win = 1; % en seconde
len_win = 2;   % en seconde

hz = 1:data_eeg.fsample/2;

quel_chan = 12;
theta_beta = [];
time_sample  = [];

for xi_win = 1 : data_eeg.fsample*slide_win : length(data_eeg.trial{1})-(len_win*data_eeg.fsample)
    
    time_sample(end+1) = xi_win;
    
    % est ce que le début de la fenetre est dans une periode rejetées
    % ou si la fin de la fenetre est dans une periode rejetées
    if sum(sum([xi_win > artifact(:,1) xi_win < artifact(:,2)],2)==2) < 1 || sum(sum([xi_win+(len_win*data_eeg.fsample) > artifact(:,1) xi_win < artifact(:,2)],2)==2) < 1
        
        ff = abs(fft(data_eeg.trial{1}(quel_chan,xi_win:xi_win+(len_win*data_eeg.fsample)), data_eeg.fsample));
        %theta_beta(end+1) = sum(ff(4:8))/sum(ff(12:40));
        theta_beta(end+1) = sum(log(ff(4:8)))/sum(log(ff(12:40)));
        
    else 
        theta_beta(end+1) = NaN;
        
    end
end, clear xi_win

output.time       = time_sample;
output.theta_beta = theta_beta;

function artifact_EOG = search_EOG_artefact(eeg_filename, trl)


% EOG
cfg     = [];
cfg.trl = [trl.sample(1),...
           trl.sample(2),...
           0];
cfg.dataset  = eeg_filename;
cfg.headerfile = eeg_filename;
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = 'EO+';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = 'yes';

[~, artifact_EOG] = ft_artifact_zvalue(cfg);

artifact_EOG = artifact_EOG - trl.sample(1);
   
function artifact_muscle = search_muscle_artefact(eeg_filename, trl) 

cfg     = [];
cfg.trl = [trl.sample(1),...
           trl.sample(2),...
           0];
cfg.dataset  = eeg_filename;
cfg.headerfile = eeg_filename;
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'F4';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [90 120];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_muscle] = ft_artifact_zvalue(cfg);

artifact_muscle = artifact_muscle - trl.sample(1);

function data = artefact_eog(data_ica, data_to_correct, chan_eog, corr_eog_th, titre)

% artefact correction based on ICA
%
% data_ica : datas no epoched, low pass-band 20 Hz
%
global GUI

% select eog datas for correlation 
cfg         = [];
cfg.channel = chan_eog.label;
data_eog = ft_selectdata(cfg, data_ica);


cfg        = [];
cfg.method = 'runica'; 
comp = ft_componentanalysis(cfg, data_ica);



if strcmp(get(findobj('Tag', titre), 'String'), 'auto')
    
    cfg = [];
    cfg.layout = 'EEG1020.lay';
    cfg.viewmode = 'component';
    ft_databrowser(cfg, comp)
    title(titre)
    
    % In Neurofeedback data, there is no EOG channel
    if ~isempty(data_eog)
        
        % correlation with EOG
        % To select the better componnent
        c = [];
        for xi = 1 : size(comp.trial{1}, 1)
            [a, b] =  corr([comp.trial{1}(xi,~isnan(comp.trial{1}(xi,:)))' data_eog.trial{1}(1,~isnan(data_eog.trial{1}))']);
            c(xi) = a(1,2);
        end, clear xi
        
        quel_comp_rej = find(abs(c) > corr_eog_th);
        
        if isempty(quel_comp_rej)
            quel_comp_rej = find(abs(c) == max(abs(c)));
        end
        
    else
    % when ther is no EOG channel, we can't choose components    
        quel_comp_rej = [];
        
    end
    
else
% if ICA selection is not auto    
    quel_comp_rej = str2num(get(findobj('Tag', titre), 'String'));
    
end

switch titre
    case 'ICA Open Eyes'
        GUI.RTB.output.ICA_OY = quel_comp_rej;
    case 'ICA Close Eyes'
        GUI.RTB.output.ICA_CY = quel_comp_rej;
end

% Correction
cfg = [];
cfg.component = quel_comp_rej;
data = ft_rejectcomponent(cfg, comp, data_to_correct);

function display_preproc(meth_artefact, th_abs_CY, th_abs_OY, win_CY, win_OY, chan_front, chan_occi, chan_eog, data_front_CY_clean, data_front_OY_clean)

global GUI 

data_front_OY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_OY, double(GUI.BLAST_Object.Fs), 0, 'no');
data_front_CY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_CY, double(GUI.BLAST_Object.Fs), 0, 'no');
temp = chan_eog; chan_eog = {}; chan_eog.label = temp;
clear temp
data_eog_CY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_CY, double(GUI.BLAST_Object.Fs), 0, 10);
data_eog_OY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_OY, double(GUI.BLAST_Object.Fs), 0, 10);

figure('Position', [100,100, 800 600], 'Color', GUI.Colors(1,:),...
     'Tag', 'data_correct')
subplot(2,1,1)
hold on
title('Close Eyes')
plot(data_eog_CY.time{1}, data_eog_CY.trial{1} - 200, 'color', [.5 .5 .5])
subplot(2,1,2)
hold on
title('Open Eyes')
plot(data_eog_OY.time{1}, data_eog_OY.trial{1} - 200, 'color', [.5 .5 .5])

for xi_chan = 1 : size(data_front_CY.trial{1},1)
    subplot(2,1,1)
    hold on
    plot(data_front_CY.time{1}, data_front_CY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
    
    subplot(2,1,2)
    hold on
    plot(data_front_CY.time{1}, data_front_OY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
end

subplot(2,1,1)
ylim([-300 200+(100*xi_chan)])
subplot(2,1,2)
ylim([-300 200+(100*xi_chan)])


if meth_artefact ~= 5
    
    win = ones(1, 1 * GUI.BLAST_Object.Fs);
    
    [xx, yy] = find(conv2(abs(data_front_CY.trial{1}) > th_abs_CY, win, 'same'));
    data_front_CY.trial{1}(xx,yy) = nan;
    [xx, yy] = find(conv2(abs(data_front_OY.trial{1}) > th_abs_CY, win, 'same'));
    data_front_OY.trial{1}(xx,yy) = nan;
    clear xx yy
    
    for xi_chan = 1 : size(data_front_CY.trial{1},1)
        subplot(2,1,1)
        hold on
        plot(data_front_CY.time{1}, data_front_CY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
        
        subplot(2,1,2)
        hold on
        plot(data_front_OY.time{1}, data_front_OY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
    end
    
    
else
    
    for xi_chan = 1 : size(data_front_CY.trial{1},1)
        subplot(2,1,1)
        hold on
        plot(data_front_CY_clean.time{1}, data_front_CY_clean.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
        
        subplot(2,1,2)
        hold on
        plot(data_front_OY_clean.time{1}, data_front_OY_clean.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
    end
       
    
end

function display_preproc_V02(win_CY, win_OY, chan_front, chan_occi, chan_eog, data_front_CY_clean, data_front_OY_clean)

global GUI 

data_front_OY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_OY, double(GUI.BLAST_Object.Fs), 0, 45);
data_front_CY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_CY, double(GUI.BLAST_Object.Fs), 0, 45);

data_eog_CY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_CY, double(GUI.BLAST_Object.Fs), 0, 10);
data_eog_OY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_OY, double(GUI.BLAST_Object.Fs), 0, 10);



figure('Position', [100,100, 1000 1000], 'Color', GUI.Colors(1,:),...
       'Tag', 'data_correct')
subplot(2,1,1)
hold on
title(['Close Eyes - channels ' strjoin(data_front_CY.label) ' + EOG'],...
      'FontSize', 17, 'Color', [.8 .8 .8])
plot(data_eog_CY.time{1}, data_eog_CY.trial{1} - 200, 'color', [.5 .5 .5])
xlim([0 data_eog_CY.time{1}(end)])
subplot(2,1,2)
hold on
title(['Open Eyes - channels ' strjoin(data_front_OY.label) ' + EOG'],...
      'FontSize', 17, 'Color', [.8 .8 .8])
plot(data_eog_OY.time{1}, data_eog_OY.trial{1} - 200, 'color', [.5 .5 .5])
xlim([0 data_eog_OY.time{1}(end)])

% plot raw datas
for xi_chan = 1 : size(data_front_CY.trial{1},1)
    subplot(2,1,1)
    hold on
    plot(data_front_CY.time{1}, data_front_CY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
    
    subplot(2,1,2)
    hold on
    plot(data_front_CY.time{1}, data_front_OY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
end

subplot(2,1,1)
ylim([-300 200+(100*xi_chan)])
subplot(2,1,2)
ylim([-300 200+(100*xi_chan)])



for xi_chan = 1 : size(data_front_CY.trial{1},1)
        subplot(2,1,1)
        hold on
        plot(data_front_CY_clean.time{1}, data_front_CY_clean.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
        
        subplot(2,1,2)
        hold on
        plot(data_front_OY_clean.time{1}, data_front_OY_clean.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
end
   
function display_preproc_V03(win_CY, win_OY, chan_front, chan_eog, data_front_CY_clean, data_front_OY_clean)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fonction plot the channels used for RTB
% After and Befor :
%           - artefact removal
%           - ICA correction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global GUI 

% Import datas  
switch GUI.Source
    
    case 'TRC_Blast'
        
        data_front_OY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_OY, double(GUI.BLAST_Object.Fs), 0, 45);
        data_front_CY = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_front, win_CY, double(GUI.BLAST_Object.Fs), 0, 45);
        
        data_eog_CY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_CY, double(GUI.BLAST_Object.Fs), 0, 10);
        data_eog_OY   = Extract_datas_Clinic_TRC(GUI.BLAST_Object.Path_TRC, chan_eog, win_OY, double(GUI.BLAST_Object.Fs), 0, 10);

    case 'MAT_Neurofeedback'
        
        data_front_OY = Extract_datas_from_neurofeedback_mat('OE', chan_front, 0, 45);
        data_front_CY = Extract_datas_from_neurofeedback_mat('CE', chan_front, 0, 45);
    
        data_eog_CY.trial{1} = {};
        data_eog_CY.label = {};
        data_eog_OY.trial{1} = {};
        data_eog_OY.label = {};

    case 'Vamp_Neurofeedback'
   
        data_front_OY = Extract_datas_from_neurofeedback_vamp(chan_front, 'open', 0, 45);
        data_front_CY = Extract_datas_from_neurofeedback_vamp(chan_front, 'close', 0, 45);

        data_eog_OY = Extract_datas_from_neurofeedback_vamp(chan_eog, 'open', 0, 10);
        data_eog_CY = Extract_datas_from_neurofeedback_vamp(chan_eog, 'close', 0, 10);

end


figure('Position', [100,100, 1000 1000], 'Color', GUI.Colors(1,:),...
       'Tag', 'data_correct')

   subplot(2,1,1)
hold on
title(['Close Eyes - channels ' strjoin(data_front_CY.label) ' + ' strjoin(data_eog_CY.label)],...
      'FontSize', 17, 'Color', [.8 .8 .8])
xlim([0 data_front_CY.time{1}(end)])

subplot(2,1,2)
hold on
title(['Open Eyes - channels ' strjoin(data_front_OY.label) ' + ' strjoin(data_eog_OY.label)],...
      'FontSize', 17, 'Color', [.8 .8 .8])
xlim([0 data_front_OY.time{1}(end)])

% plot raw datas
for xi_chan = 1 : size(data_front_CY.trial{1},1)
    subplot(2,1,1)
    hold on
    plot(data_front_CY.time{1}, data_front_CY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
    
    subplot(2,1,2)
    hold on
    plot(data_front_OY.time{1}, data_front_OY.trial{1}(xi_chan, :) + ((xi_chan - 1)*200), 'r')
    
end

subplot(2,1,1)
if ~isempty(data_eog_CY.trial{1})
    plot(data_eog_CY.time{1}, data_eog_CY.trial{1} - 300, 'color', [.5 .5 .5])
    ylim([-300-(abs(nanmin(data_eog_CY.trial{1}(:)))) 300*xi_chan])
else
    ylim([-300 300*xi_chan])

end
subplot(2,1,2)
if ~isempty(data_eog_OY.trial{1})
    plot(data_eog_OY.time{1}, data_eog_OY.trial{1} - 300, 'color', [.5 .5 .5])
    ylim([-300-(abs(nanmin(data_eog_OY.trial{1}(:)))) 300*xi_chan])
else
    ylim([-300 300*xi_chan])
end
clear xi_chan



cfg         = [];
cfg.channel = chan_front.label;
data_front_OY_clean = ft_selectdata(cfg, data_front_OY_clean);
data_front_CY_clean = ft_selectdata(cfg, data_front_CY_clean);


% plot trial keeped
for xi_trial = 1 : numel(data_front_CY_clean.trial)
    for xi_chan = 1 : size(data_front_CY_clean.trial{xi_trial},1)
        subplot(2,1,1)
        hold on
        plot(data_front_CY_clean.time{xi_trial} + data_front_CY_clean.sampleinfo(xi_trial,1)/data_front_CY.fsample, data_front_CY_clean.trial{xi_trial}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
    end, clear xi_chan
end, clear xi_trial

for xi_trial = 1 : numel(data_front_OY_clean.trial)
    for xi_chan = 1 : size(data_front_OY_clean.trial{xi_trial},1)
        subplot(2,1,2)
        hold on
        plot(data_front_OY_clean.time{xi_trial} + data_front_OY_clean.sampleinfo(xi_trial,1)/data_front_OY.fsample, data_front_OY_clean.trial{xi_trial}(xi_chan, :) + ((xi_chan - 1)*200), 'k')
    end, clear xi_chan
end, clear xi_trial

function Display_output()

global GUI


% Display spectrum mean and frequency bands
freq_id = nearest(GUI.RTB.output.ff.CY.mean.freq, 1): nearest(GUI.RTB.output.ff.CY.mean.freq, 25);
fft_diff = squeeze(nanmean(GUI.RTB.output.ff.CY.mean.powspctrm(:,freq_id),1)) - squeeze(nanmean(GUI.RTB.output.ff.OY.mean.powspctrm(:,freq_id),1));
figure('Position', [100,100, 500 500],  'MenuBar', 'no', 'Color', GUI.Colors(1,:),...
       'Tag', 'Frequency_power');
fill([GUI.RTB.output.theta_band GUI.RTB.output.theta_band(end:-1:1)],...
     [min(fft_diff) min(fft_diff) max(fft_diff*1.3) max(fft_diff*1.3)], 'k', 'FaceAlpha', .25, 'EdgeAlpha', .3)
hold on
fill([GUI.RTB.output.beta_band GUI.RTB.output.beta_band(end:-1:1)],...
     [min(fft_diff) min(fft_diff) max(fft_diff*1.3) max(fft_diff*1.3)], 'k', 'FaceAlpha', .25, 'EdgeAlpha', .3)
plot(GUI.RTB.output.ff.CY.mean.freq(freq_id), squeeze(nanmean(GUI.RTB.output.ff.CY.mean.powspctrm(:,freq_id),1)), '.--b', 'LineWidth', 2);
plot(GUI.RTB.output.ff.OY.mean.freq(freq_id), squeeze(nanmean(GUI.RTB.output.ff.OY.mean.powspctrm(:,freq_id),1)), '.--k', 'LineWidth', 2)
plot(GUI.RTB.output.ff.CY.mean.freq(freq_id), fft_diff(freq_id), 'r', 'LineWidth', 2)
plot([GUI.RTB.output.alpha_peak GUI.RTB.output.alpha_peak], [min(fft_diff) max(fft_diff*1.3)], '-k', 'LineWidth', 3)
legend({'Theta', 'Beta', 'Eyes Close', 'Eyes Open', 'Difference', 'Alpha Pic'})
xlabel('Hz')
clear fft_diff




RTB_OY = GUI.RTB.output.RTB.trial.OY;
RTB_CY = GUI.RTB.output.RTB.trial.CY;


mean_cum_CY = RTB_CY(1);
for xi_trial = 2 : length(RTB_CY)
   mean_cum_CY = [mean_cum_CY median(RTB_CY(1:xi_trial))];
end
mean_cum_OY = RTB_OY(1);
for xi_trial = 2 : length(RTB_OY)
   mean_cum_OY = [mean_cum_OY median(RTB_OY(1:xi_trial))];
end

% Display
max_RTB = max([mean(RTB_CY), mean(RTB_OY)]) + 3*max([std(RTB_CY), std(RTB_OY)]);
min_RTB = min([mean(RTB_CY), mean(RTB_OY)]) - max([std(RTB_CY), std(RTB_OY)]);

figure('Position', [100,100, 800 900],  'MenuBar', 'no', 'Color', GUI.Colors(1,:),...
       'Tag', 'RTB_time'),
subplot(3,2,1), plot(RTB_OY, 'k')
hold on, fill([1, length(RTB_OY), length(RTB_OY), 1],...
               [nanmedian(RTB_OY)+nanstd(RTB_OY), nanmedian(RTB_OY)+nanstd(RTB_OY), nanmedian(RTB_OY)-nanstd(RTB_OY), nanmedian(RTB_OY)-nanstd(RTB_OY)],...
            'y', 'FaceAlpha', 0.2) 
          plot([1 length(RTB_OY)],...
               [nanmedian(RTB_OY) nanmedian(RTB_OY)], 'Color', [.2 .6 1], 'LineWidth', 2)
          plot([1 length(RTB_OY)],...
               [GUI.RTB.output.RTB.mean.OY GUI.RTB.output.RTB.mean.OY], 'g', 'LineWidth', 2, 'Color', [0 .6 .3])
          text(3, 0.9*max_RTB, 'RTB', 'FontSize', 12)
          text(10, 0.8*max_RTB, ['by windows : ', num2str(nanmedian(RTB_OY)), ' (-/+ ', num2str(nanstd(RTB_OY)), ')'],...
               'FontSize', 12, 'Color', [.2 .6 1])
          text(10, 0.7*max_RTB, ['spect mean : ', num2str(GUI.RTB.output.RTB.mean.OY)],...
               'FontSize', 12, 'Color', [0 .6 .3])
           title('RTB - Yeux Ouverts', 'FontSize', 15, 'Color', [.8 .8 .8])
           xlabel('time (window)')
           ylim([min_RTB, max_RTB])
           xlim([0 length(RTB_OY)])
           
subplot(3,2,2), plot(RTB_CY, 'g')
hold on, fill([1, length(RTB_CY), length(RTB_CY), 1],...
               [nanmedian(RTB_CY)+nanstd(RTB_CY), nanmedian(RTB_CY)+nanstd(RTB_CY), nanmedian(RTB_CY)-nanstd(RTB_CY), nanmedian(RTB_CY)-nanstd(RTB_CY)],...
            'y', 'FaceAlpha', 0.2) 
          plot([1 length(RTB_CY)],...
               [nanmedian(RTB_CY) nanmedian(RTB_CY)], 'Color', [.2 .6 1], 'LineWidth', 2)
          plot([1 length(RTB_CY)],...
               [GUI.RTB.output.RTB.mean.CY GUI.RTB.output.RTB.mean.CY], 'g', 'LineWidth', 2, 'Color', [0 .6 .3])
          text(3, 0.9*max_RTB, 'RTB', 'FontSize', 12)
          text(10, 0.8*max_RTB, ['by windows : ', num2str(nanmedian(RTB_CY)), ' (-/+ ', num2str(nanstd(RTB_CY)), ')'],...
               'FontSize', 12, 'Color', [.2 .6 1])
          text(10, 0.7*max_RTB, ['spect mean : ', num2str(GUI.RTB.output.RTB.mean.CY)],...
               'FontSize', 12, 'Color', [0 .6 .3])
           title('RTB - Yeux Fermés', 'FontSize', 15, 'Color', [.8 .8 .8])
           xlabel('time (window)')
           ylim([min_RTB, max_RTB])
           xlim([0 length(RTB_CY)])
clear max_RTB min_RTB



max_RTB = max([mean(mean_cum_CY), mean(mean_cum_OY)]) + 2*max([std(mean_cum_CY), std(mean_cum_OY)]);
min_RTB = min([mean(mean_cum_CY), mean(mean_cum_OY)]) - max([std(mean_cum_CY), std(mean_cum_OY)]);

subplot(3,2,3), plot(mean_cum_OY, 'k')
hold on
plot([1 length(mean_cum_OY)], [mean_cum_OY(end) mean_cum_OY(end)], 'k', 'LineWidth', 2)
plot([1 length(RTB_OY)],...
     [GUI.RTB.output.RTB.mean.OY GUI.RTB.output.RTB.mean.OY], 'g', 'LineWidth', 2, 'Color', [0 .6 .3])
ylim([min_RTB, max_RTB])
title('Mediane cumulee - Yeux Ouverts', 'FontSize', 15, 'Color', [.8 .8 .8])
xlim([0 length(RTB_OY)])

subplot(3,2,4), plot(mean_cum_CY, 'k')
hold on
plot([1 length(mean_cum_CY)], [mean_cum_CY(end) mean_cum_CY(end)], 'k', 'LineWidth', 2)
plot([1 length(RTB_CY)],...
     [GUI.RTB.output.RTB.mean.CY GUI.RTB.output.RTB.mean.CY], 'g', 'LineWidth', 2, 'Color', [0 .6 .3])
ylim([min_RTB, max_RTB])
title('Mediane cumulee - Yeux Fermes','FontSize', 15, 'Color', [.8 .8 .8])
xlim([0 length(RTB_CY)])




max_alpha = max([mean(GUI.RTB.output.alpha.trial.CY), mean(GUI.RTB.output.alpha.trial.OY)]) + 2*max([std(GUI.RTB.output.alpha.trial.CY), std(GUI.RTB.output.alpha.trial.OY)]);
min_alpha = min([mean(GUI.RTB.output.alpha.trial.CY), mean(GUI.RTB.output.alpha.trial.OY)]) - max([std(GUI.RTB.output.alpha.trial.CY), std(GUI.RTB.output.alpha.trial.OY)]);

subplot(3,2,5), plot(GUI.RTB.output.alpha.trial.OY, 'k')
ylim([min_alpha, max_alpha])
title('Alpha', 'FontSize', 15, 'Color', [.8 .8 .8])
xlim([0 length(RTB_OY)])

subplot(3,2,6), plot(GUI.RTB.output.alpha.trial.CY, 'k')
ylim([min_alpha, max_alpha])
title('Alpha', 'FontSize', 15, 'Color', [.8 .8 .8])
xlim([0 length(RTB_CY)])



% Screan
fprintf('\n\n\nFeatures:\n\n')
fprintf('Channels used to compute RTB:\n');
for xi_chan = 1 : numel(GUI.RTB.param_artefact.chan_front.String)
    fprintf('\t\t\t\t%s\n', GUI.RTB.param_artefact.chan_front.String{xi_chan});
end, clear xi_chan

fprintf('\nThreshold\n')
fprintf('\tYeux ouverts : %s microVolt\n ', num2str(GUI.RTB.output.threshold_OY))
fprintf('\tYeux fermes  : %s microVolt\n ', num2str(GUI.RTB.output.threshold_CY))

fprintf('\nICA composantes\n')
fprintf('\tYeux ouverts : %s \n ', num2str(GUI.RTB.output.ICA_OY))
fprintf('\tYeux fermes  : %s \n ', num2str(GUI.RTB.output.ICA_CY))

fprintf('\n\nResultats:\n')
fprintf('Pourcentage de signal utile\n')
fprintf('\tYeux ouverts : %s\n ', num2str(round(GUI.RTB.output.pct_OY*100)))
fprintf('\tYeux fermes  : %s\n ', num2str(round(GUI.RTB.output.pct_CY*100)))

fprintf('\nBande de frequence\n')
fprintf('\talpha : %s Hz\n ', num2str(GUI.RTB.output.alpha_peak))
fprintf('\ttheta : [%s %s] Hz\n ', num2str(GUI.RTB.output.theta_band(1)),...
                                   num2str(GUI.RTB.output.theta_band(2)))
fprintf('\tbeta  : [%s %s] Hz\n ', num2str(GUI.RTB.output.beta_band(1)),...
                                   num2str(GUI.RTB.output.beta_band(2)))

fprintf('\nRatio theta/beta\n')
fprintf('Par fenetre; \t       mediane (deviation standard)\n')
fprintf('\tYeux ouverts : %s (+/- %s)\n', num2str(median(GUI.RTB.output.RTB.trial.OY)),...
                                           num2str(std(GUI.RTB.output.RTB.trial.OY)))
fprintf('\tYeux fermes  : %s (+/- %s)\n', num2str(median(GUI.RTB.output.RTB.trial.CY)),...
                                           num2str(std(GUI.RTB.output.RTB.trial.CY)))
fprintf('Moyenne des fenetres; \n')
fprintf('\tYeux ouverts : %s \n ', num2str(GUI.RTB.output.RTB.mean.OY))
fprintf('\tYeux fermes  : %s \n ', num2str(GUI.RTB.output.RTB.mean.CY))
                                  
function Save_output()

global GUI

% Select the output directory
[dir_raw,~,~] = fileparts(GUI.BLAST_Object.Path_TRC);
GUI.RTB.output.dir_output = uigetdir(dir_raw, 'Pick a Directory to OUTPUT');
clear dir_raw

[~,filename,~] = fileparts(GUI.BLAST_Object.Path_TRC);

% save output .mat file
output = GUI.RTB.output;
save(fullfile(GUI.RTB.output.dir_output, [filename '_metrics.mat']),...
     'output')
clear output


% Display
p_signal = findobj('Tag', 'data_correct');
p_RTB = findobj('Tag', 'RTB_time');

saveas(p_signal, fullfile(GUI.RTB.output.dir_output, [filename '_Signal']), 'jpg')
saveas(p_signal, fullfile(GUI.RTB.output.dir_output, [filename '_Signal']), 'fig')
saveas(p_RTB, fullfile(GUI.RTB.output.dir_output, [filename '_RTB']), 'jpg')
saveas(p_RTB, fullfile(GUI.RTB.output.dir_output, [filename '_RTB']), 'fig')

fid = fopen(fullfile(GUI.RTB.output.dir_output, [filename '.txt']), 'w+');

fprintf(fid, '%s\n\n', filename);
% Screan
fprintf(fid, 'Features:\n\n');
fprintf(fid, 'Channels used to compute RTB:\n');
for xi_chan = 1 : numel(GUI.RTB.param_artefact.chan_front.String)
    fprintf(fid, '\t\t\t\t%s\n', GUI.RTB.param_artefact.chan_front.String{xi_chan});
end, clear xi_chan
fprintf(fid, '\nThreshold\n');
fprintf(fid, '\tYeux ouverts : %s microVolt\n ', num2str(GUI.RTB.output.threshold_OY));
fprintf(fid, '\tYeux fermes  : %s microVolt\n ', num2str(GUI.RTB.output.threshold_CY));

fprintf(fid, '\nICA composantes\n');
fprintf(fid, '\tYeux ouverts : %s \n ', num2str(GUI.RTB.output.ICA_OY));
fprintf(fid, '\tYeux fermes  : %s \n ', num2str(GUI.RTB.output.ICA_CY));

fprintf(fid, '\n\nResultats:\n');
fprintf(fid, 'Pourcentage de signal utile\n');
fprintf(fid, '\tYeux ouverts : %s\n ', num2str(round(GUI.RTB.output.pct_OY*100)));
fprintf(fid, '\tYeux fermes  : %s\n ', num2str(round(GUI.RTB.output.pct_CY*100)));

fprintf(fid, '\nBande de frequence\n');
fprintf(fid, '\ttheta : [%s %s] Hz\n ', num2str(GUI.RTB.output.theta_band(1)),...
                                   num2str(GUI.RTB.output.theta_band(2)));
fprintf(fid, '\tbeta  : [%s %s] Hz\n ', num2str(GUI.RTB.output.beta_band(1)),...
                                   num2str(GUI.RTB.output.beta_band(2)));

fprintf(fid, '\nRatio theta/beta\n')
fprintf(fid, 'Par fenetre; \t       mediane (deviation standard)\n')
fprintf(fid, '\tYeux ouverts : %s (+/- %s)\n', num2str(median(GUI.RTB.output.RTB.trial.OY)),...
                                           num2str(std(GUI.RTB.output.RTB.trial.OY)))
fprintf(fid, '\tYeux fermes  : %s (+/- %s)\n', num2str(median(GUI.RTB.output.RTB.trial.CY)),...
                                           num2str(std(GUI.RTB.output.RTB.trial.CY)))
fprintf(fid, 'Moyenne des fenetres; \n')
fprintf(fid, '\tYeux ouverts : %s \n ', num2str(GUI.RTB.output.RTB.mean.OY))
fprintf(fid, '\tYeux fermes  : %s \n ', num2str(GUI.RTB.output.RTB.mean.CY))

fclose(fid);
