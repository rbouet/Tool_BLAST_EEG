function [HDR, EVT, DATA] = ahdr2mat(ahdrfilename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [HDR, EVT, DATA] = ahdr2mat(ahdrfilename)
%
% This function import EEG file 
%           V-amp (*.ahdr)
%
% Author: MB
%
% Add FEALDTRIP path
%
% Author : MM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath(genpath('C:\_MANU\_U821\_DirMatlab\fieldtrip'))

[Fpath,Fname,~] = fileparts(ahdrfilename);
TMP_vhdrfilename   = [Fpath,filesep,Fname,'.vhdr'];
TMP_vmrkfilename   = [Fpath,filesep,Fname,'.vmrk'];
eegfilename        = [Fpath,filesep,Fname,'.eeg'];

% Copy *.eeg and *.amrk files
copyfile([Fpath,filesep,Fname,'.amrk'], TMP_vmrkfilename);
copyfile([Fpath,filesep,Fname,'.ahdr'], TMP_vhdrfilename);

% Build Header
HDR = ft_read_header(TMP_vhdrfilename);
HDR_tmp = HDR;
HDR_tmp.nChans=HDR.nChans + 1;
HDR_tmp.orig.NumberOfChannels = HDR_tmp.nChans;
HDR_tmp.nSamples =  HDR.nSamples*HDR.nChans/(HDR.nChans + 1);
HDR_tmp.orig.nSamples =HDR_tmp.nSamples;

% Import datas
begsample = 1;
endsample = HDR_tmp.nSamples;
chanindx = 1:HDR_tmp.nChans-1;
data = ft_read_data(eegfilename,'header',HDR_tmp,'begsample',begsample,'endsample',endsample,'chanindx',1:HDR.nChans);

DATA = data';
ix2cut = unique(sum(~DATA(end-100:end,:)));
DATA(end-ix2cut+1:end,:) = [];

HDR.nSamples = size(DATA,1);
HDR.orig.nSamples = size(DATA,1);


% Import events
[event] = ft_read_event(TMP_vmrkfilename);
sampevt = cell2mat({event.sample});

itmp = 1 ;
for ievt = 1:length(event)
    if (strfind(event(ievt).type,'New Segment'))
        ixnewseg(itmp) = ievt;
        itmp = itmp + 1;
    end
end

if (~isempty(char(event.value)))
    sampevt(ixnewseg)=[];
    valevt = char(event.value);
    valevt(ixnewseg,:)=[];
    valevt_tmp = valevt;
    valevt_tmp(:,1)=[];
    
    if isempty(str2num(valevt_tmp))
        EVT = {sampevt',valevt};
    else
        EVT = [sampevt',str2num(valevt_tmp)];
    end
else
    EVT = [];
end


delete(TMP_vmrkfilename)
delete(TMP_vhdrfilename)