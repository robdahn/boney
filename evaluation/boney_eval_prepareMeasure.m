

    % plot bone measure for different sites
   % out.bone_med_tiv{1} = out.bone_med{1} %./ (out.TIV{1}/median(out.TIV{1}(:)))'; normalise for TIV
    site = 8;
 
    if all(isnan(out.rGMV{1})), paravol = {}; else, paravol = {'rGMV';'rWMV';'rCMV';'TIV'}; end
    para =  [{
      'bmdh';
      'bmdt';
      'bmdf'; 
      }
    paravol; ...
      {
      'age'; 
      'VAT';'ASAT';'fatp';'bmi';'waist'; 
      'age_meno';'HRT';
      'sittingheight'; ...
      'birthweight';... 'sunburnchild';
      'timesummer';'walk';'modPA';'vigPA'; ...
      'happyness';'healthsat';... 'timewinter';
      'alc';'smoking';'fIQ' ... 'smokingn';
      }];
    para2 =  {
      'bmdh';'bmdt';'bmdf';'age';...
      'VAT';'ASAT';'fatp';'bmi';'waist';'birthweight';'sittingheight';...
      'smoking';'timesummer';'happyness';'healthsat';...'smokingN';'timewinter';
      'walk';'modPA';'vigPA';'alc'; ...
      'sunburnchild'};

out.BMDhead{1} = bmdh'; out.BMDtotal{1} = bmdt'; out.BMDfemur{1} = bmdf'; 
    


switch matonly
  case 0
    measure = { 
      'bone_num'; 
      'bone_med'; 
      'bone_std'; 
      'spm_bone_med1'; 
      'spm_bone_med2'; 
      'spm_bone_med3';
    }; 
  case 1 
    measure = { 
      'bone_med'; 
      'spm_bone_med1'; 
      'spm_bone_med2'; 
      'spm_bone_med3';
    }; 
  case 4
    measure = { 
      'bone_med';
    }; 
  case 5
    measure = { 
      'bone_num'; 
      'bone_med'; 
      'bone_std'; 
    };
  case 10
    measure = { 
      'bone_num'; 
      'bone_med'; 
    };
  case 11 %{11,12}
    measure = { 
      'tis_bonecortex';
      'tis_bonemarrow'; ... 'tis_bonedensity';
      ... 'bone_num';   ... bone volume 
      'bone_med';       ... classical bone med measure of tissue intensity in MRI ('bone_std' is similar)
      'tis_bone';       ... weighted average SPM head intensity peaks
      'tis_head';       ... weighted average SPM head intensity peaks
      'tismri_volfatr'; ... volume with high intensity (fatvolume)
      ... 'tismri_volmusr'; % volume of low intensy head structures (weak) 
      'vROI_bonecortex3'; 'vROI_bonethickness3'; 'vROI_headthickness3'; % simple short version
      ... 'vROI_bonecortex2'; 'vROI_bonecortex3'; 'vROI_bonecortex4'; 
      ... 'vROI_bonemarrow2'; 'vROI_bonemarrow3'; 'vROI_bonemarrow4';
      ... 'vROI_bonethickness2'; 'vROI_bonethickness3'; 'vROI_bonethickness4'; 
      ... 'vROI_headthickness2'; 'vROI_headthickness3'; 'vROI_headthickness4'; 
      }; 
    if matonly == 12
      measure = [measure ; {
        'sROI_bonecortex3'; 'sROI_bonethickness3'; 'sROI_headthickness3';
        ... 'sROI_bonecortex2'; 'sROI_bonecortex3'; 'sROI_bonemarrow4'; 
        ... 'sROI_bonemarrow2'; 'sROI_bonemarrow3'; 'sROI_bonemarrow4';
        ... 'sROI_bonethickness2'; 'sROI_bonethickness3'; 'sROI_bonethickness4'; 
        ... 'sROI_headthickness2'; 'sROI_headthickness3'; 'sROI_headthickness4'; 
      }];
    end 
  case 12
    switch plevel
      case 0
        % == Super short only with final measures == 
        % Thickness measures are defined on the segments and therefore less 
        % dependent on the image properties but intensities support more 
        % subtil grading. 
        measurex = { 
          ... 1) Bone measures
          'tis_bone'              0   % weighted average of the SPM bone intensity peaks
          'sROI_bonecortex3'      1   % final bone-cortex intensity measure
          'sROI_bonethickness3'   0   % final bone thickness measure
          'sROI_ibonethickness3'  1   % final bone thickness measure
          'sROI_BMDH'             0   % final combined measure
           ...
          ... 2) Head measures
          'tis_head'              2   % weighted average SPM head intensity peaks
          'tismri_volfatr'        1   % final relative fat-volume (fat/total head volume)
          'sROI_headthickness3'   0   % final head thickness measure 
          ...
          'BMDhead'               2
          'BMDtotal'              0
          'BMDfemur'              0
          ...
          'rGMV'                  2          
          'rWMV'                  0
          'rCMV'                  0
        }; 
      case 1
        % a bit longer with historical and parameter steps
        measurex = { 
          ... 1) Bone measures 
          'tis_bone'              0   % weighted average of the SPM bone intensity peaks (important!)
          ... 1.1) Extract own value - used median to reduce outliers 
     ...     'bone_med'              1   % classical bone med measure of tissue intensity
          ... 1.2) Surface-based extraction on the global bone measure to focus on surface measures 
          %'vROI_bone1'            1    
          %'sROI_bone1'            0   % >>   
          ... 1.3) Local separation of bone cortex and bone marrow 
     ...     'sROI_bonecortex1'      1   % >>
     ...     'sROI_bonemarrow1'      0
          ... 1.4) Region dependencies for major lobes to focus on occiptial regions 
          'sROI_bonecortex1'      1
          'sROI_bonecortex2'      0
          'sROI_bonecortex3'      0  % >>
          'sROI_bonecortex4'      0
          ... 1.5) Bone thickness
          'sROI_ibonecortex3'     1
          %'sROI_bonethickness1'   1  
          %'sROI_bonethickness2'   0  
          'sROI_bonethickness3'   0  
          %'sROI_bonethickness4'   0  
          'sROI_BMDH'             1   % final combined measure
          ...
          ... 2) Head measures
          'tis_head'              2   % weighted average SPM head intensity peaks
          'tismri_volfatr'        1   % final relative fat-volume (fat/total head volume)
          'sROI_headthickness3'   0   % final head thickness measure 
          ...'shdt3'                 0   % final combined head measure
          ...
          ... 3) Combined measures ?
          ...'sROI_nbhthickness3'    2
          ...
          'BMDhead'               2
          'BMDtotal'              0
          'BMDfemur'              0
          ...
          'rGMV'                  2          
          'rWMV'                  0
          'rCMV'                  0
         }; 
      case 2
        measurex = { 
          %'tis_bonecortex'        0   % minimum SPM bone intensity peak
          %'tis_bonemarrow'        0   % maximum SPM bone intensity peak
          'tis_bone'              0   % weighted average of the SPM bone intensity peaks
          %'tis_head'                 1 % weighted average of the SPM head intensity peaks
          %'tis_bonecortexpercentage' 0 % SPM Gaussian weights - low correlation
          %'tis_headfatpercentage'    0 % SPM Gaussian weights - low correlation 
          ...
          %'bone_med'              2   % classical bone med measure of tissue intensity in MRI ('bone_std' is similar)
          ...
          %'vROI_bone1'            2
          %'vROI_bone2'            0
          %'vROI_bone3'            0 
          %'vROI_bone4'            0
          %'sROI_bone1'            1
          %'sROI_bone2'            0
          %'sROI_bone3'            0 
          %'sROI_bone4'            0
          ...
          'vROI_bonecortex1'      2
          'vROI_bonecortex2'      0
          'vROI_bonecortex3'      0
          'vROI_bonecortex4'      0
          'sROI_bonecortex1'      1
          'sROI_bonecortex2'      0
          'sROI_bonecortex3'      0 
          'sROI_bonecortex4'      0
          'sROI_ibonecortex1'      1
          'sROI_ibonecortex2'      0
          'sROI_ibonecortex3'      0 
          'sROI_ibonecortex4'      0
          ...
          %'vROI_bonemarrow1'      2
          %'vROI_bonemarrow2'      0
          %'vROI_bonemarrow3'      0
          %'vROI_bonemarrow4'      0
          %'sROI_bonemarrow1'      1
          %'sROI_bonemarrow2'      0
          %'sROI_bonemarrow3'      0 
          %'sROI_bonemarrow4'      0
          ...
          'vROI_bonethickness1'   2  % simple short version
          'vROI_bonethickness2'   0 
          'vROI_bonethickness3'   0  
          'vROI_bonethickness4'   0  
          'sROI_bonethickness1'   1  
          'sROI_bonethickness2'   0  
          'sROI_bonethickness3'   0  
          'sROI_bonethickness4'   0  
          %'sROI_nbonethickness1'   1  % this add a TIV bias
          %'sROI_nbonethickness2'   0  
          %'sROI_nbonethickness3'   0  
          %'sROI_nbonethickness4'   0  
          ...
          'vROI_headthickness1'   2
          'vROI_headthickness2'   0
          'vROI_headthickness3'   0
          'vROI_headthickness4'   0
          'sROI_headthickness1'   1
          'sROI_headthickness2'   0
          'sROI_headthickness3'   0
          'sROI_headthickness4'   0
          %'sROI_nheadthickness1'  1
          %'sROI_nheadthickness2'  0
          %'sROI_nheadthickness3'  0
          %'sROI_nheadthickness4'  0
          ...
          %'sROI_bhthickness1'     2
          %'sROI_bhthickness2'     0
          %'sROI_bhthickness3'     0
          %'sROI_bhthickness4'     0
          %'sROI_nbhthickness1'    1
          %'sROI_nbhthickness2'    0
          %'sROI_nbhthickness3'    0
          %'sROI_nbhthickness4'    0
          ...
          %'vROI_BMDH'             2   % combined measures
          %'sROI_BMDH'             0   % combined measures
          %'vROI_BMDH1'            1   % combined measures
          %'vROI_BMDH2'            0   % combined measures
          %'vROI_BMDH3'            0   % combined measures
          %'vROI_BMDH4'            0   % combined measures
          ...
          %'sROI_nbonethickness3'  2  
          %'sROI_headthickness3'   1
          %'sROI_nheadthickness3'  0
          %'sROI_nbhthickness3'    1
          ...
          'tis_head'              2   % weighted average SPM head intensity peaks
          'tismri_volfatr'        1   % volume with high intensity (fatvolume)
          %'vhdt1'                 0   % = tis_head ???
          %'shdt2'                 0
          %'shdt3'                 0
          ...
          %'vROI_bonemix'          2
          %'vROI_boneM2C'          0
          %'vROI_boneMxT'          0
          %'sROI_bonemix'          1       
          %'sROI_boneM2C'          0
          %'sROI_boneMxT'          0
          ...
          'BMDhead'               2
          'BMDtotal'              0
          'BMDfemur'              0
          ...
          'rGMV'                  2          
          'rWMV'                  0
          'rCMV'                  0
          }; 
    end
    measure = measurex(:,1);
    nyi  = find(cell2mat(measurex(:,2)')==1)-1; 
    nyim = find(cell2mat(measurex(:,2)')==2)-1; % [ 6 14 22 28 36 numel(measure)-6 numel(measure)-3 ]; % short

  case {3,2}
    %out.iBMDhead{1} = -bmdh'; out.iBMDtotal{1} = -bmdt'; out.iBMDbody{1} = -bmdf'; 
    measure  = { ...
  ... bone intensities (surface and volume based are quite similar)
      ... globeal bone intensities (mean ~ median, std ~ iqr) >  
      1 'bone_med';       %  classical measure (korrelation to TIV that makes it better?)
      1 'boneskmeans21';  %  = 'mnBone' ... better than 31
      1 'boneskmeans22';  %  = 'mnMarrow' this is best for rFAT 
      1 'boneskmeans33';  %  = this is best for ASAT (esp. in femals)
  ... thickness
      1 'bone_vol';
      1 'hdt';  
      1 'b2ht'; 
  ... volumes 
      1 'bonehdt2';  
      1 'bonehdt6'; 
  ...
      1 'BMDhead'; 
      1 'BMDtotal';
      1 'BMDbody';
      1 'rGMV'; 
      1 'rWMV'; 
      1 'rCMV';
    }; 
end
%measure = measurep(measurep(:,1)>plevel,2);