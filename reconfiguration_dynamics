%======================================================================
%          APPLY RECONFIGURATION DYNAMICS WITH SPECIFIC ROIS
%======================================================================
%@author: Marloes Bet
%@email:  m.bet@amsterdamumc.nl
%last updated: 25 11 2024 
%status: FINISHED
%to do: Let Tommy review code

%Review History
%NONE YET

% Description:
% Determines community assignment across dynamic windows based on inputted timeseries. Outputs flexibility,
% promiscuity, cohesion and disjointedness per node and averaged across all nodes corresponding to a specific monoaminergic system.
% Assumes that scripts for calculating flexibility, promiscuity, and cohesion matrix are available in working directory (otherwise add to path).
%
% Run:
% - fMRI preprocessing
% - steps to determine which region belongs to which network
% (literature-based or based on available data)
% 
% Input:
% - fMRI timeseries
% - general assignment of regions to cortical networks, and to monoaminergic systems
% 
% Output:
% - reconfiguration dynamics variables for each node separately as well as
% averaged across monoaminergic systems
%-----------------------------------------------------------------------

%% ESTABLISH COMMUNITIES
addpath('/scratch/anw/mdbet/fatigue/scripts')

% Get original RSN assignment
tmp=zeros(224,1);
% Default Mode Network (DMN):
tmp([3,4,5,6,11,12,13,14,23,27,33,35,39,41,42,43,44,51,52,79,81,82,...
    83,84,87,88,95,113,121,122,123,141,143,144,151,153,154,165,175,...
    176,178,179,181,182,187,188])=1;
% Frontoparietal Network (FPN):
tmp([16,17,18,19,20,21,22,24,26,28,29,30,31,32,34,36,40,99,100,137,...
    138,142,147,177])=2;
% Dorsal Attention Network (DAN):
tmp([7,8,25,55,56,63,64,85,86,91,92,97,98,107,125,126,127,128,129,...
    130,133,134,136,139,140,148,150])=3;
% Ventral Attention Network (VAN):
tmp([1,2,15,37,38,61,62,124,145,146,149,166,167,168,169,170,173,174,...
   180,183,184,185,186])=4;
% Visual Network (VN):
tmp([105,106,108,114,119,120,135,152,189,190,191,192,193,194,195,196,...
    197,198,199,200,201,202,203,204,205,206,207,208,209,210])=5;
% Somatomotor Network (SMN):
tmp([9,10,53,54,57,58,59,60,65,66,67,68,71,72,73,74,75,76,80,131,132,...
    155,156,157,158,159,160,161,162,163,164,171,172])=6;
% Limbic System (LS): (usually excluded due to distortions in fMRI data)
tmp([45,46,47,48,49,50,69,70,77,78,89,90,93,94,96,101,102,103,104,109,...
    110,111,112,115,116,117,118])=7;
% Deep Gray Matter (DGM) / subcortical:
tmp([211,212,213,214,215,216,217,218,219,220,221,222,223,224])=8; %Subcortical
% Cerebellum (excluded due to distortions in fMRI data)
%tmp(225)=9;

% Community assignment based on 5% highest receptor density
HT1a=[110,109,113,114,117,215,222,103,111,104,94,169];
HT2a=[81,85,87,97,88,82,99,39,37,123,86,23];
HTT=[213,220,211,218,214,221];
HT_all=unique([HT1a,HT2a,HTT]);
DA=[220,213,219,212,224,217];
DAT=[217,220,214,224,213,221,212,219];
DA_all=unique([DA,DAT]);
NAT=[211,218,67,68,162,59,57,72,155,161,156,60];

%% DYNAMIC COMMUNITY ASSIGNMENT
% Load timeseries per atlas region
timeseries=dir('/scratch/anw/mdbet/fatigue/timeseries_Graz/*_BNA_timeseries.txt');
timeseries=timeseries(~ismember({timeseries.name},{'.','..'}));

% Initiate output file
output=cell(length(timeseries),5);
count=0;

if exist('scratch/anw/mdbet/fatigue/output_data/community_assignments.mat') == 2
    disp('Community assignment information has already been calculated.');
    load('/scratch/anw/mdbet/fatigue/output_data/community_assignments.mat')
else
    disp('Now calculating community assignment for each dynamic window...');

    for ii = timeseries' %loop over rows in the struct
        count=count+1;
    
        % Match subject IDs to the participant list in the output file
        subjectid = split(ii.name,'_');
        subjectid = subjectid{1};
    
        disp("Processing part "+count+": "+subjectid)
    
        % Place subject id in output matrix
        output{count,1}=subjectid;
        
        % Calculate functional connectivity matrix
        ts=dlmread([ii.folder filesep ii.name]);
        fMRIconmatrix=abs(atanh(corr(ts)));%correlate timeseries, Fishers r-to-z transformation, and make absolute;
        fMRIconmatrix(logical(eye(size(fMRIconmatrix))))=NaN; %convert diagonal line values to 0; EC formula cannot handle NaN or Inf values
        output{count,2}=fMRIconmatrix;
    
        % Determine communities in individual subjects' dynamic time windows based on Yeo cortical network assignment.
        output{count,3}=community_detection(ts,tmp,27,1);
    
    end
    save('community_assignments.mat', 'output');
end

%% CALCULATING RECONFIGURATION VARIABLES FOR ALL REGIONS
% Load the cell matrix with community assignment if it is not yet present
% in the current workspace.
if exist('output', 'var') == 1
    disp('The cell matrix "output" exists in the workspace. Moving on with this file...');
else
    disp('The cell matrix "output" does not exist in the workspace. Loading from files...');
    load('/scratch/anw/mdbet/fatigue/output_data/community_assignments.mat','output')
end

% Loop over the cell matrix's rows to calculate variables for each subject.
count=0;

for ii = 1:length(output)
    pp=output{ii,1};
    disp(pp);

    assignment=output{ii,3};
    S=assignment';

    % Calculate node flexibility and promiscuity
    output{ii,4}=flexibility(S);
    output{ii,5}=promiscuity(S);

    % Calculate cohesion and disjointedness and related variables
    [cohesion_mat,node_cohesion,node_disjointedness,node_flexibility,cohesion_strength]=calc_node_cohesion(assignment);

%    Cij: Cohesion matrix, similar to adjacency matrix but edge
%                weights represent the number of times two nodes change
%                communities mutually
    output{ii,6}=cohesion_mat;

%    Node_cohesion: Node cohesion for a given node
%            
%                             # Times Node Changes Communities Mutually
%             Cohesiveness = -------------------------------------------
%                             # Possible Times Nodes Change Communities
    output{ii,7}=node_cohesion;

%    Node_disjoint: Node disjointedness for a given node
%            
%                            # Times Node Changes Communities Independently
%           Disjointedness = ----------------------------------------------
%                              # Possible Times Nodes Change Communities
    output{ii,8}=node_disjointedness;

%    Node_flexibility: Node flexibility for a given node
%            
%                                # Times Node Changes Communities
%              Flexibility = -------------------------------------------
%                             # Possible Times Nodes Change Communities
    output{ii,9}=node_flexibility; % should be the same matrix as in output{ii,4}


end

save('reconfiguration_output.mat','output')

%% AVERAGING VARIABLES ACROSS REGIONS OF INTEREST PER MONOAMINERGIC RECEPTOR OR TRANSPORTER
if exist('output', 'var') == 1
    disp('The cell matrix "output" exists in the workspace. Moving on with this file...');
else
    disp('The cell matrix "output" does not exist in the workspace. Loading from files...');
    load('/scratch/anw/mdbet/fatigue/reconfiguration_output.mat','output');
end

for ii=1:length(output)
    pp=output{ii,1};
    disp(pp);

    % Access node vectors for each variable
    flex_vec=output{ii,4};
    prom_vec=output{ii,5};
    cohe_vec=output{ii,7};
    disj_vec=output{ii,8};

    % Average flexibility per monoaminergic receptor/transporter system
    %output{ii,10}=mean(flex_vec(D1));
    output{ii,11}=mean(flex_vec(DA));
    output{ii,12}=mean(flex_vec(DAT));
    output{ii,13}=mean(flex_vec(DA_all));
    output{ii,14}=mean(flex_vec(HT1a));
    output{ii,15}=mean(flex_vec(HT2a));
    output{ii,16}=mean(flex_vec(HTT));
    output{ii,17}=mean(flex_vec(HT_all));
    output{ii,18}=mean(flex_vec(NAT));
    
    % Average promiscuity per monoaminergic receptor/transporter system
    %output{ii,19}=mean(prom_vec(D1));
    output{ii,20}=mean(prom_vec(DA));
    output{ii,21}=mean(prom_vec(DAT));
    output{ii,22}=mean(prom_vec(DA_all));
    output{ii,23}=mean(prom_vec(HT1a));
    output{ii,24}=mean(prom_vec(HT2a));
    output{ii,25}=mean(prom_vec(HTT));
    output{ii,26}=mean(prom_vec(HT_all));
    output{ii,27}=mean(prom_vec(NAT));

    % Average cohesion per monoaminergic receptor/transporter system
    %output{ii,28}=mean(cohe_vec(D1));
    output{ii,29}=mean(cohe_vec(DA));
    output{ii,30}=mean(cohe_vec(DAT));
    output{ii,31}=mean(cohe_vec(DA_all));
    output{ii,32}=mean(cohe_vec(HT1a));
    output{ii,33}=mean(cohe_vec(HT2a));
    output{ii,34}=mean(cohe_vec(HTT));
    output{ii,35}=mean(cohe_vec(HT_all));
    output{ii,36}=mean(cohe_vec(NAT));

    % Average disjointedness per monoaminergic receptor/transporter system
    %output{ii,37}=mean(disj_vec(D1));
    output{ii,38}=mean(disj_vec(DA));
    output{ii,39}=mean(disj_vec(DAT));
    output{ii,40}=mean(disj_vec(DA_all));
    output{ii,41}=mean(disj_vec(HT1a));
    output{ii,42}=mean(disj_vec(HT2a));
    output{ii,43}=mean(disj_vec(HTT));
    output{ii,44}=mean(disj_vec(HT_all));
    output{ii,45}=mean(disj_vec(NAT));

    % Average variables across whole brain
    output{ii,46}=mean(flex_vec,'all');
    output{ii,47}=mean(prom_vec,'all');
    output{ii,48}=mean(cohe_vec,'all');
    output{ii,49}=mean(disj_vec,'all');
    
    % Correct values for whole-brain value per variable
    for jj=11:18
        % 10 t/m 18 or 50 t/m 58: flexibility
        output{ii,jj+40}=output{ii,jj}/output{ii,46}; 
    end

    for jj=20:27
        % 19 t/m 27 or 59 t/m 67: promiscuity
        output{ii,jj+40}=output{ii,jj}/output{ii,47}; 
    end

    for jj=29:36
        % 28 t/m3 6 or 68 t/m 76: cohesion
        output{ii,jj+40}=output{ii,jj}/output{ii,48}; 
    end

    for jj=38:45
        % 37 t/m 45 or 77 t/m 85: disjointedness
        output{ii,jj+40}=output{ii,jj}/output{ii,49}; 
    end


end

% Create column names for each value 10-85.
colNames = {'ID','D1_flex','D1D2_flex','DAT_flex','DA_flex','5HT1a_flex','5HT2a_flex','5HTT_flex','5HT_flex','NAT_flex',...
    'D1_prom','D1D2_prom','DAT_prom','DA_prom','5HT1a_prom','5HT2a_prom','5HTT_prom','5HT_prom','NAT_prom',...
    'D1_cohe','D1D2_cohe','DAT_cohe','DA_cohe','5HT1a_cohe','5HT2a_cohe','5HTT_cohe','5HT_cohe','NAT_cohe',...
    'D1_disj','D1D2_disj','DAT_disj','DA_disj','5HT1a_disj','5HT2a_disj','5HTT_disj','5HT_disj','NAT_disj',...
    'global_flex','global_prom','global_cohe','global_disj',...
    'D1_flex_corr','D1D2_flex_corr','DAT_flex_corr','DA_flex_corr','5HT1a_flex_corr','5HT2a_flex_corr','5HTT_flex_corr','5HT_flex_corr','NAT_flex_corr',...
    'D1_prom_corr','D1D2_prom_corr','DAT_prom_corr','DA_prom_corr','5HT1a_prom_corr','5HT2a_prom_corr','5HTT_prom_corr','5HT_prom_corr','NAT_prom_corr',...
    'D1_cohe_corr','D1D2_cohe_corr','DAT_cohe_corr','DA_cohe_corr','5HT1a_cohe_corr','5HT2a_cohe_corr','5HTT_cohe_corr','5HT_cohe_corr','NAT_cohe_corr',...
    'D1_disj_corr','D1D2_disj_corr','DAT_disj_corr','DA_disj_corr','5HT1a_disj_corr','5HT2a_disj_corr','5HTT_disj_corr','5HT_disj_corr','NAT_disj_corr'};

% Add column names to the top of the data matrix
outputWithHeaders = [colNames; output(:,[1,10:end])];

% Save variables of interest as Excel file
writetable(cell2table(outputWithHeaders(2:end,:), 'VariableNames', outputWithHeaders(1,:)), '/scratch/anw/mdbet/fatigue/output_data/reconf_vars_per_receptor_incl_corrected.xlsx', 'WriteRowNames', true);


