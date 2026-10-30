function out = boney_segment_simpleBone(seg8t, Ym, Yp, refine, vx_vol)
%simpleBone. First simple version that worked surprisingly well.
% There is a basic correction of the SPM bone segmentation to avoid 
% misclassification of high-intensity bone marrow as head tissue. However, 
% there is no further masking of critical regions (e.g. facial bones) or 
% evaluation of different regions.
% 
%  out = boney_segment_simpleBone(seg8t,Ym,Yp,refine,vx_vol)
%
%  seg8t      .. SPM unified segmentation structure with tissue peaks etc.
%  Ym         .. intensity normalized image
%  Yp         .. SPM tissue classes
%  refine     .. do additional refinement
%  vx_vol     .. voxel size (of the resolution reduced images)
%  out        .. output structure with basic information
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
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


  if ~exist('vx_vol','var'), vx_vol = [1 1 1]; end
  vxcm3 = prod(vx_vol) / 1000; % voxel volume in cm3

  % load SPM tissue segments (the BlaiottaTPM classes are already mapped)
  for ci = 1:(max(seg8t.lkp) - 1) 
     out.int_ci(1,ci) = median(Ym(Yp{ci}>0.9));       % intensity
     out.vol_ci(1,ci) = sum(Yp{ci}(:)>0.5) * vxcm3;   % volume
  end
  % estimate TIV
  out.TIV = sum(out.vol_ci(1,1:3));


  % create bone mask
  % As the SPM segmentation often misaligns high values of the bone marrow  
  % to the head, a correction by morphological operations is required. 
  if refine
    Ymsk = smooth3(Yp{4})>0.25;                   % remove small noise
    Ymsk = cat_vol_morph(Ymsk,'l');               % get the largest bones element
    Ymsk = cat_vol_morph(Ymsk,'dc',6);            % remove holes in the bone
  else
    Ymsk = Yp{4}>0.5;  
  end  

  % extract CSF-normalized values (hist)
  out.con       = out.int_ci(1,3);                  % CSF contrast
  out.bone_num  = (sum(Ymsk(:)) * vxcm3) / out.TIV; % bone volume normalised by TIV
  out.bone_med  = median(Ym(Ymsk(:)))   / out.con;  % median bone intensity value
  out.bone_std  = std(Ym(Ymsk(:)))      / out.con;  % standard deviation of the bone intensity                  
end
