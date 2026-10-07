%% Config file
%{

Written by Victoria St Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

Set your study-specific paths & variables to run the xGLM pipeline. 

%}

% Paths
analysisPath = 'your\path\here';

% All reg labels
reglabels = {'collab'; 'ind';'scaffolding';
    'partnerLeftPFC_collab';'partnerLeftPFC_ind';
    'partnerRightPFC_collab';'partnerRightPFC_ind';
    'partnerLeftTPJ_collab';'partnerLeftTPJ_ind';
    'partnerRightTPJ_collab';'partnerRightTPJ_ind';'constant'};

% Behavioural regs
beh_regs = {'cns_handing';'cns_positioning';'cns_pointing';'cns_demonstrating';
    'cns_correcting';'cws_handing';'cws_positioning';'cws_pointing';
    'cws_demonstrating'; 'cws_correcting'
};

% number ROIs
numROIs = 4;

% roi Names
roiNames = {'leftPFC';'rightPFC';'leftTPJ';'rightTPJ'};

% Channels belonging to each region of interest
leftPFC = [1 2 4 5];
rightPFC = [6 7 8 9];
leftTPJ = [10 11 12 13];
rightTPJ = [14 15 17 18];

% Channel numbers corresponding to short-separation channels
sscs = [3; 16];

% NIRS sampling rate
fs = 25;
