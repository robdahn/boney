% Bone aging test script.
% * processing of different datasets such as IXI, CAMCAN and NKI datasets 
%   of healthy (adults) subjects with full age range
%
% * Evaluation of the current version by (not including image data but CSV files?) ... 
%   * BMD measures in the UKB (non public by UKB)
%   * Hammers data?
%   * IXI / NKI / CamCAN >> aging effects (public)
%   * OASIS >> CT vs. MRI , longitudinal, test-retest (public)
%   * Buchert, Simon >> large-scale test-retest  (public)
%   * small test subdataset (UKB360)
%
% * Evaluation of selected parameters of multiple version 
%   * selected small dataset?
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


%
% ToDo: 
%  * restructure these batches and add parts of the data and description 
%    how to use (eg. download IXI and set path) 
%  * include and organize different test conditions and subdirs with
%    detailed reports 
%  * how to include the version or date?
%  * extend help 
%  * include setting for our data paths

clear groups

% data - groups (bias corrected m*.nii)
% * fasttest includes a smaller two group test subset of young and old IXI
%   subjects.
% * matonly==1 uses only the SPM preprocessing parameter of the *_seg8.mat 
%   (named as SPM) whereas matonly==0 call the so called MED pipeline that 
%   loads and refines the SPM tissue segments to estimate own volume and
%   intensity parameters.
% 0-MED, 1-SPMo;   2-SPMn, 3-MED2, 4-SPMo0, 5-MED0;    10-final0, 11-final1, 12-final2 
matonlyn  = {'MED','SPM',   'NEWSPM','NEWMED','SPM0','MED0','','','','','BoneySPM','BoneyV','BoneyS'}; 
if 1
  maindirs  = {
    '/Users/robertdahnke/MRData/Boney/'; 
    '/Users/robertdahnke/MRData/Boney/'; 
    '/Users/robertdahnke/MRData/Boney/'; 
    };
elseif true
  maindirs  = {
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    '/Volumes/SG5TB/MRData/202301_MRIbones/'; 
    };
else
  maindirs  = {
    '/Users/polona/Downloads/202301_MRIbones';
    '/Volumes/OneTouch/MRData/202301_MRIbones';
    '/Users/polona/Downloads/IXI'; 
    }; 
end
resultdir = fullfile(maindirs{1},'results_Boney_20230828'); 
if ~exist(resultdir,'dir'), mkdir(resultdir); end
dname = 'aging';

% get the image data
%testset = % IXI, UKB
groups{1} = [ 
    ...cat_vol_findfiles(maindirs{3},'mIXI*.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'camcan','BIDSsep','anat','msub*T1w.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'NKI_RAW','NKI_T1','msub*T1w.nii');
    ...NA - cat_vol_findfiles(maindirs{2},'ukb_bones','msub*T1w.nii');
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_testset')      ,'mG*.nii'  ,struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_testsetExt')   ,'mG*.nii'  ,struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100a'),'msub*.nii',struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100b'),'msub*.nii',struct('maxdepth',1)));
    cellstr(cat_vol_findfiles(fullfile(maindirs{1},'ukb_bones_subsample100c'),'msub*.nii',struct('maxdepth',1)));
    %cellstr(cat_vol_findfiles(fullfile(maindirs{1},'OASIS3_CT')    ,'c04_*CTseg.nii',struct('maxdepth',4)));
    ... setdiff( ...
    ...    cat_vol_findfiles(maindirs{2},'OASIS3','BIDS','msub*T1w.nii'), ...
    ...    cat_vol_findfiles(maindirs{2},'OASIS3','BIDS','msub*echo*T1w.nii'));
    ...setdiff( ...
      ...cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','m*T1w.nii'), ...
      ...[cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','mw*T1w.nii');
       ...cat_vol_findfiles(maindirs{2},'NKI_RS','T1w','ME*T1w.nii')]);
    ]; 

% get the phenotypic data
Pcsv{1} = {
    ...fullfile(maindirs{3},'IXI.csv');
    ...fullfile(maindirs{2},'camcan','BIDSsep','cc700-scored','participant_data.csv'); 
    ...fullfile(maindirs{2},'NKI_RAW','participants.tsv'); 
    ...fullfile(maindirs{2},'NKI_RAW',NKI_RS/201506/NKI.1-39.phenotypic.csv';
    %fullfile(maindirs{1},'OASIS3_CT','+tables/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv');
    fullfile(maindirs{1},'ukb_bones_testset','ukb_3000.csv');
    };

%
boney_eval_readCSV

cat_get_defaults('extopts.expertgui',2)


% single tests
if 0
%% Processing
  matonly = 12;   
  out1 = boney_getBoneMeasures({groups{1}(1)},struct('matonly',matonly,'rerun',0));
  %out = getBoneMeasures(groups,struct('matonly',matonly,'rerun',0));
end

%%
%for matonly = [12, 11, 10, 0:5] % 0-MED, 1-SPMo;   2-SPMn, 3-MED2, 4-SPMo0, 5-MED0;    10-final0, 11-final1, 12-final2 
matonly = 12;   
out = boney_getBoneMeasures(groups,struct('matonly',matonly,'rerun',0));
% end
%%
plevel = 2; % 1-basic, 2-3-full
boney_eval_prepareMeasure

% scatter plot for each test paramter with linear fits
% print overview    
printfg = 0; 
boney_eval_printDetails % (this also involve some other steps >> seperate!)
boney_eval_printMain

%% print for each (bone) measure of one bone processing pipeline a detailed
printfg = 1;
boney_eval_printDetails
   
%%
boney_eval_printCrossCorrelation



%% create table
tab = cell(1,numel(Pcsv{1}));
for csvi = 1:numel(Pcsv{1})
   tab{csvi}{1,1}       = 'Filename'; 
   tab{csvi}(2:numel(groups{csvi})+1,1) = spm_str_manip(groups{csvi},'t');

   [pp,ff,ee] = spm_fileparts(groups{1}{si}); 
   if contains(pp,'OASIS3')
       subj = groups{csvi}; 
       
       tab{csvi}(1,2:4) = {'OASISID','DAY','RUN'}; 
       for si = 1:numel(subj)
           tab{csvi}{si+1,2} = spm_str_manip(subj{si},'hhht'); tab{csvi}{si+1,2} = tab{csvi}{si+1,2}(5:end);
           tab{csvi}{si+1,3} = spm_str_manip(subj{si},'hht');  tab{csvi}{si+1,3} = tab{csvi}{si+1,3}(6:end); 
           if contains(spm_str_manip(subj{si},'t'),'run')
             tab{csvi}{si+1,4} = spm_str_manip(subj{si},'t'); tab{csvi}{si+1,4} = tab{csvi}{si+1,4}(strfind(tab{csvi}{si+1,4},'run-')+[4 5]); 
           else 
             tab{csvi}{si+1,4} = '01';
           end
       end
   end
   
   
   FN = {'bone_num','bone_med','rGMV','rWMV','rCMV'}; 
   for fni = 1:numel(FN)
       tab{csvi}{1,end+1} = FN{fni};
       if size( out.(FN{fni}){csvi} , 2) > 1  
          tab{csvi}(2:numel(groups{csvi})+1,end) = num2cell(out.(FN{fni}){csvi}');
       else
          tab{csvi}(2:numel(groups{csvi})+1,end) = num2cell(out.(FN{fni}){csvi});
       end
   end
   
   [pp,ff,ee] = fileparts(Pcsv{1}{csvi}); 
   resultcsv{csvi} = fullfile(pp,sprintf('%s_bonemeasures.csv',ff));
   cat_io_csv(resultcsv{csvi},tab{csvi});
end


