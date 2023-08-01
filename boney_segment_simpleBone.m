function out = boney_segment_simpleBone(seg8t,Ym,Yp)
%simpleBone. First simple version that worked quite well.

% == (1) SPM gaussians ==
  % get high prob value of CSF tissue
  out.spm_con = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
    ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
  
  if max(seg8t.lkp) > 6
    out.spm_bone_med1 = seg8t.mn( seg8t.lkp == 4 ) / out.spm_con;
    out.spm_bone_med2 = seg8t.mn( seg8t.lkp == 5 ) / out.spm_con;
  else
    % default case for 6 TPM classes
    bone = sort(seg8t.mn( seg8t.lkp == 4 ) / out.spm_con); % order intensities within class 4 (normalized by CSF)
    out.spm_bone_med1 = bone(min(numel(bone),3)); % highest value
    out.spm_bone_med2 = bone(1); 
    out.spm_bone_med3 = bone(min(numel(bone),2)); 
  end
 

% == (2) median bone intensity ==
  % load SPM tissue segments
  for ci = 1:5 
     out.int_ci(1,ci) = median(Ym(Yp{ci}>0.9));     % intensity
     out.vol_ci(1,ci) = sum(Yp{ci}(:)>0.5) / 1000;  % volume
  end
  % estimate TIV
  out.TIV = sum(out.vol_ci(1,1:3));

  % create bone mask
  if max(seg8t.lkp) > 6
    Ymsk = smooth3(Yp{4} + Yp{5})>0.5;          % remove small noise
  else
    Ymsk = smooth3(Yp{4})>0.25;                 % remove small noise
  end
  Ymsk = cat_vol_morph(Ymsk,'l');               % get the largest bones element
  Ymsk = cat_vol_morph(Ymsk,'dc',6);            % remove holes in the bone
  
  % extract values (hist)
  out.con       = out.int_ci(1,3);                  % CSF contrast
  out.bone_num  = (sum(Ymsk(:)) / 1000) / out.TIV;  % bone volume normalised by TIV
  out.bone_med  = median(Ym(Ymsk(:))) / out.con;    % median bone intensity value
  out.bone_std  = std(Ym(Ymsk(:)))    / out.con;                     
end