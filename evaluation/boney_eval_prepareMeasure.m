

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
      'birthweight';'sunburnchild';
      'timesummer';'timewinter';'happyness';'healthsat';...
      'smoking';'fIQ' ...
      }];
    para2 =  {
      'bmdh';'bmdt';'bmdf';'age';...
      'VAT';'ASAT';'fatp';'bmi';'waist';'birthweight';'sittingheight';...
      'smoking';'timesummer';'timewinter';'happyness';'healthsat';...
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
      case 1
        measure = { 
          'sROI_bonecortex3'; ... regionbased new 
          'vROI_bonemarrow3';    
          'vROI_BMDH';              % combined measures
          ...
          'vhdt1'; 
          'sROI_headthickness3';
          ...
          'BMDhead'; 
          'BMDtotal';
          'BMDbody';
          'rGMV'; 
          'rWMV'; 
          'rCMV'; 
          }; 
          nyi  = [ 3 ]; 
          nyim = [ numel(measure)-6 numel(measure)-3 ]; % short
      case 2
        measure = { 
          'tis_bonecortex'; 
          'tis_bonemarrow'; ... 'tis_bonedensity';
          'tis_bone';       ... weighted average SPM head intensity peaks
          'bone_med';       ... classical bone med measure of tissue intensity in MRI ('bone_std' is similar)
          'vROI_bonecortex3';    
          'sROI_bonecortex3';
          ...
          'vROI_bonemarrow0';    
          'vROI_bonemarrow1';    
          'vROI_bonemarrow2';    
          'vROI_bonemarrow3';    
          'vROI_bonemarrow4';    
          'sROI_bonemarrow1'; 
          'sROI_bonemarrow2'; 
          'sROI_bonemarrow3'; 
          'sROI_bonemarrow4'; 
          ...
          'vROI_bonethickness1';    % simple short version
          'sROI_bonethickness1';    
          'vROI_bonethickness2';    % simple short version
          'sROI_bonethickness2';    
          'vROI_bonethickness3';    % simple short version
          'sROI_bonethickness3';    
          'vROI_bonethickness4';    % simple short version
          'sROI_bonethickness4';    
          ...
          'vROI_BMDH';              % combined measures
          'vROI_BMDH2';              % combined measures
          'vROI_BMDH1';              % combined measures
          'vROI_BMDH3';              % combined measures
          'vROI_BMDH4';              % combined measures
          'sROI_BMDH'; 
          'tis_head';       ... weighted average SPM head intensity peaks
          'tismri_volfatr'; ... volume with high intensity (fatvolume)
          ...
          'sROI_headthickness1';  
          'sROI_headthickness2';  
          'sROI_headthickness3';
          'sROI_headthickness4';
          'vhdt1'; 
          'shdt3'; 
          ...
          'BMDhead'; 
          'BMDtotal';
          'BMDfemur';
          'rGMV'; 
          'rWMV'; 
          'rCMV'; 
          }; 
          nyi  = [ 3 6 ]; 
          nyim = [ 7 15 23 29 numel(measure)-6 numel(measure)-3 ]; % short
    end
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