function out = boney_segment_simpleBone(seg8t,Ym,Yp)
%simpleBone. First simple version that worked surprisingly good.
% There is a basic correction of the SPM bone segmentation to avoid miss-
% classification of high-intensity bone marrow as head tissue. However, 
% there is no further masking of critical regions (e.g. face bones) or 
% evaluation of different regions.
% 
%  out = boney_segment_simpleBone(seg8t,Ym,Yp)
%
%  seg8t      .. SPM unified segmentation structure with tissue peaks etc.
%  Ym         .. intensity normalized image
%  Yp         .. SPM tissue classes
%  out        .. output structure with basic information
%   .spm_con  .. 
%   .spm_bone_med# .. bone values just given by SPM tissue peaks 
%   .int_ci   .. median tissue intensity of each tissue class
%   .vol_ci   .. volume of each tissue class
%   .TIV      .. total intracranial volume
%   .con      .. CSF contrast used for bone intensity normalization 
%   .bone_num .. relative bone volume normalized by TIV 
%   .bone_med .. median bone intensity 
%   .bone_std .. standard deviation of bone intensity 
%
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
% $Id$



  % == (1) SPM gaussians per tissue class ==
  % get high prob value of CSF tissue
  out.spm_con = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
    ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
  
  if max(seg8t.lkp) > 6
    % RD202308: two-class bone tissue
    % This was a special test case with artificial seperation of the bone
    % tissue class into inner bonemarrow and outer hard bone. However, I am
    % not sure if the classes are correctly used in following functions,  
    % i.e. this case should maybe removed!
    cat_io_cprintf('warn','\nWarning: a-typical TPM test-case that should have a soft and a hard bone class!\n')
    out.spm_bone_med1 = seg8t.mn( seg8t.lkp == 4 ) / out.spm_con;
    out.spm_bone_med2 = seg8t.mn( seg8t.lkp == 5 ) / out.spm_con;
  else
    % default case for 6 TPM classes and typical/various SPM subtissues 
    % (i.e., [ 1 1 2 3 4 2 ]) with the typical 3 peaks for bone
    bone = sort(seg8t.mn( seg8t.lkp == 4 ) / out.spm_con); % order intensities within class 4 (normalized by CSF)
    out.spm_bone_med1 = bone(min(numel(bone),3)); % highest value (i.e. bone marrow) 
    out.spm_bone_med2 = bone(1);                  % lowest values (i.e. hard bone) 
    out.spm_bone_med3 = bone(min(numel(bone),2)); % second median value if available otherwise highest value
    clear bone
  end
 


  
  % == (2) median bone intensity ==
  % load SPM tissue segments
  for ci = 1:(max(seg8t.lkp) - 1) 
     out.int_ci(1,ci) = median(Ym(Yp{ci}>0.9));     % intensity
     out.vol_ci(1,ci) = sum(Yp{ci}(:)>0.5) / 1000;  % volume
  end
  % estimate TIV
  out.TIV = sum(out.vol_ci(1,1:3));


  % create bone mask
  % As the SPM segmentation often miss-align high values of the bone marrow  
  % to the head, a correction by morphological options is required. 
  if max(seg8t.lkp) > 6
    Ymsk = smooth3(Yp{4} + Yp{5})>0.5;          % remove small noise
  else
    Ymsk = smooth3(Yp{4})>0.25;                 % remove small noise
  end
  Ymsk = cat_vol_morph(Ymsk,'l');               % get the largest bones element
  Ymsk = cat_vol_morph(Ymsk,'dc',6);            % remove holes in the bone
  

  % extract CSF-normalized values (hist)
  out.con       = out.int_ci(1,3);                  % CSF contrast
  out.bone_num  = (sum(Ymsk(:)) / 1000) / out.TIV;  % bone volume normalised by TIV
  out.bone_med  = median(Ym(Ymsk(:)))   / out.con;  % median bone intensity value
  out.bone_std  = std(Ym(Ymsk(:)))      / out.con;  % standard deviation of the bone intensity                  
end