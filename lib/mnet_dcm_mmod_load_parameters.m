function [mC,Aw,Cw,nobs,ns,nm] = mnet_dcm_mmod_load_parameters(nodedef,netdef,U)
%[mC,Cw,Aw,nobs]=set_parameters(nodedef,netdef)
% input: 
% nodedef=node.csv
% netdef=connect.csv
% U: external input
%
% output:
% mC  : model parameters
% Cw  : initial C (non-zero value is require to consider the parameter)
% Aw  : initial A (non-zero value is require to consider the parameter)
% nobs: number of observation -> nbobs.vsdi, nobs.cai, nobs.bold
% ns  : number of source (total number of neural population)
% nm  : number of modularity input
%
% Jiyoung Kang: Apil-08-2019
% updated Jiyoung Kang: 2019-09-13
% last updated Jiyoung Kang: 2020-04-06

node = readnode_csv(nodedef);

% specify connectivity
% ====================
if(strcmp(netdef,'full'))
    numnode = length(node.ID);
    Aw = ones(numnode); 
    for ii = 1:numnode
        Aw(ii,ii) = 0;
    end
else
    Aw = getAw(netdef,node.ID); 
end
A = Aw ~= 0;
%Aw(Aw == 1.0) = 0.1;   % inter-region
%Aw(Aw == 1.2) = 0.15;  % intra-region


%% Convert to output type
ID = node.ID;
CwInput = node.CwInput;
Vrest = node.Vrest;
Tau = node.Tau;
Vth = node.Vth;
MaxFreq = node.MaxFreq;
La = node.La;
type1 = node.type1;
H = node.H;
%Name = node.Name;
RegionID_VSDI = node.RegionID_VSDI;
RegionID_CaI = node.RegionID_CaI;
RegionID_BOLD = node.RegionID_BOLD;
obs_VSDI = node.obs_VSDI;
obs_CaI = node.obs_CaI;
obs_BOLD = node.obs_BOLD;

Cw = CwInput;  % input from external stimulus
C = Cw ~= 0;

mC.A = A;
mC.C = C;

ns = length(ID);
nobs.vsdi = max(RegionID_VSDI);
nobs.cai = max(RegionID_CaI);
nobs.bold = max(RegionID_BOLD);



mC.H    = H;       % PSP
mC.Fr   = MaxFreq; % maximal firing rates (Hz) %mod!
mC.Vm   = Vrest;   % resting membrane potentials (mV)
mC.thr  = Vth;     % threshold membrane potential (mV)

% ExIn: 1 = excitatory, -1 = inhibitory
mC.ExIn     = zeros(ns,1);
id          = find(strcmp(type1,'E'));
mC.ExIn(id) = 1;
id          = find(strcmp(type1,'I'));
mC.ExIn(id) = -1;   

mC.Popu.vsdi = RegionID_VSDI; %[v] population indicator
mC.La        = La;            %[-] lambda for VSDI

mC.Popu.cai  = RegionID_CaI;  %[v] population indicator for CaI
mC.Popu.bold = RegionID_BOLD; %[v] population indicator for BOLD

% added
mC.obs.vsdi = obs_VSDI; %[v] obseved signal indicator 1=on, 0=off
mC.obs.cai  = obs_CaI;  %[v] obseved signal indicator
mC.obs.bold = obs_BOLD; %[v] obseved signal indicator

mC.Tau       = Tau;


% time-constant for decaying firing rate
% --------------------------------------
mC.T    = mC.Tau * 0.001; %[v] inversion of Tau (ms-1) %% NOTE: consider a sampling rate (Hz) of measuring. 
mC.R    = 2/3;            %[v] 2/3 (Moran et al., 2013, frontiers Comp. Neurosci.)
mC.A0   = 0.17;           % CaI KSJung

% For CaI    
mC.thrCa = -27.89;      % [B] [fixed] [mV] Half-activation voltage for HVA, Lee (2013)
mC.RCa  = 0.2;  % [fixed] Slope of the sigmoid activation, Ca, Rahmati et al. (2016)

mC.TCa  = 0.8;  % [unfixed] Time constant of the kernel of calcium [ver1 1/1.25]
mC.KCa  = 0.26; % [unfixed] Converting and scaling parameter % mod. [!--ver1.K2=15]

mC.GCa  = 5;    % [fixed] Maximal conductance for calcium ions, Rahmati et al. (2016)
mC.KF   = 9.85; % [B] [fixed] Scale parameter for the ﬂuorescence trace, Yaksi and Friedrich, 2006
mC.dF   = -3.2833;    %  calibrated to 0 for initial Ca signals

mC.Ca   = 100;  % [fixed] Baseline calcium ion concentration, Rahmati et al. (2016)     
mC.Kd   = 200;  % [fixed] Dissociation coefficient,  Rahmati et al. (2016)

mC.ECa  = 120;  % Reverse potential of Ca ions  



[ns,~]  = size(Aw);     % number of neural populations
nm      = size(U.u,2);  % number of modularity input

end

function [Aw] = getAw(filename, ID)
[from, to, w] = readconnect(filename);
nNode = size(ID, 1);
Aw = zeros(nNode);

for i = 1:length(from)
    fromID = find(strcmp(ID, from{i}));
    toID = find(strcmp(ID, to{i}));
    Aw(toID, fromID) = w(i);
end

end

function [from,to,w] = readconnect(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [FROM,TO,W] = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   [FROM,TO,W] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [from,to,w] = importfile('connectKnown.csv',2, 21);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/01/13 16:36:38

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
from = dataArray{:, 1};
to = dataArray{:, 2};
w = dataArray{:, 3};
end


function tbl = readnode_csv(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [ID, CWINPUT, VREST, TAU, VTH, MAXFREQ, LA, TYPE1, H, NAME,
%  REGIONID_VSDI, REGIONID_CAI, REGIONID_BOLD] = IMPORTFILE(FILENAME)
%  reads data from text file FILENAME for the default selection.
%  Returns the data as column vectors.
%
%  [ID, CWINPUT, VREST, TAU, VTH, MAXFREQ, LA, TYPE1, H, NAME,
%  REGIONID_VSDI, REGIONID_CAI, REGIONID_BOLD] = IMPORTFILE(FILE,
%  DATALINES) reads data for the specified row interval(s) of text file
%  FILENAME. Specify DATALINES as a positive scalar integer or a N-by-2
%  array of positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [ID, CwInput, Vrest, Tau, Vth, MaxFreq, La, type1, H, Name, RegionID_VSDI, RegionID_CaI, RegionID_BOLD] = importfile("/Users/jiyoungkang/calc/DCM/simulElecPhys/05vsdiCaIbold/vsdi_cai_bold_ver2p0/node.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 12-Sep-2019 16:45:43

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 16);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ID", "CwInput", "Vrest", "Tau", "Vth", "MaxFreq", "La", "type1", "H", "Name", "RegionID_VSDI", "RegionID_CaI", "RegionID_BOLD","obs_VSDI", "obs_CaI", "obs_BOLD"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "string", "double", "string", "double", "double", "double","double", "double", "double"];
opts = setvaropts(opts, [1, 8, 10], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 8, 10], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
% ID = tbl.ID;
% CwInput = tbl.CwInput;
% Vrest = tbl.Vrest;
% Tau = tbl.Tau;
% Vth = tbl.Vth;
% MaxFreq = tbl.MaxFreq;
% La = tbl.La;
% type1 = tbl.type1;
% H = tbl.H;
% Name = tbl.Name;
% RegionID_VSDI = tbl.RegionID_VSDI;
% RegionID_CaI = tbl.RegionID_CaI;
% RegionID_BOLD = tbl.RegionID_BOLD;
% obs_VSDI = tbl.obs_VSDI;
% obs_CaI = tbl.obs_CaI;
% obs_BOLD = tbl.obs_BOLD;

end