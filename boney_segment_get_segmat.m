function [ seg8t, tis, vx_vol ] = boney_segment_get_segmat(out,verb)   
%get_segmat. Get and evaluate SPM preprocessing structure.
%
%  [seg8t, tis, vx_vol] = get_spm8(P)
%
%   P       .. one image
%   seg8t   .. spm8 structure without larger fields
%   tis     .. measures based on SPM tissue thresholds and 
%              basic information based on image and tissue properties
%   vx_vol  .. voxel-size
% 
% Comments: 
% * SPM tissue peaks for [GM WM CSF bone head background]
%   - default is [1 1 2 3 4 2]: 
%     . in T1 the lower CSF peak is the right one whereas the other one is 
%       PVE to GM or in the best case meninges?  
%       > would need a separate class?
%     . the lower CSF threshold leads to overestimation in younger subjects - PVE?
%       > the CSF value is less robust for normalization 
%     . in mixed tissues, higher mg result often in higher vr that bias the mn  
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


  % TODO: 
  %  (1) lkp tests
  %  (2) skull-stripping/defacing tests
  %  (3) TIV estimate  


  % load SPM/CTseg or  mat 
  if strcmp(out.P.seg8(end-2:end),'mat')
    seg8t          = load(out.P.seg8); 
    if out.CTseg
      seg8t.dat.model.gmm = rmfield(seg8t.dat.model.gmm,{'T','Sig'});
      % #### this is not fully working and the classes values are strange ...
      seg8t.lkp     = seg8t.sett.gmm.mg_ix;
      seg8t.mn      = seg8t.dat.model.gmm.m; 
      seg8t.mg      = seg8t.dat.model.gmm.gam'; 
      seg8t.vr      = seg8t.dat.model.gmm.W; 
       
      tmp           = spm_load_priors8(fullfile(spm('dir'),'TPM','TPM.nii'));
      seg8t.tpmA    = tmp;
      seg8t.tpm     = rmfield(tmp.V,'private');
      seg8t.Affine  = eye(4);
      seg8t.image   = spm_vol(out.P.org);
      seg8t.isCTseg = 1; 
      seg8t         = rmfield(seg8t,'dat');
      seg8t.sett    = rmfield(seg8t.sett,'B');
    else
      seg8t         = rmfield(seg8t,{'Twarp','Tbias','MT'}); % store this one for later 
      seg8t.isCTseg = 0; 
      
      if numel( seg8t.image ) > 1
        iid         = contains(out.P.org,{seg8t.image.fname}); 
        seg8t.image = seg8t.image(iid);
        seg8t.mn    = seg8t.mn(iid,:);
        seg8t.vr    = seg8t.vr(iid,iid,:);
      end
    end
    vx_vol          = sqrt(sum(seg8t.image(1).mat(1:3,1:3).^2));

% ###################
% It would be helpful to include quality measures here (or later for speed up) 
% by using the CAT QC batch
% ###################
  else
    % read CAT XML rather than the SPM file 
    Sxml  = cat_io_xml(out.P.seg8); 
    seg8t = Sxml.SPMpreprocessing;
    tmp           = spm_load_priors8(fullfile(spm('dir'),'TPM','TPM.nii'));
    seg8t.tpmA    = tmp;
    seg8t.tpm     = rmfield(tmp.V,'private');
    if ~isfield(seg8t,'mg') % older cat versions don't have this field!
      for i = 1:max(seg8t.lkp)
        seg8t.mg(seg8t.lkp == i,1) = 1 / sum(seg8t.lkp == i); 
      end
    end
    seg8t.image   = spm_vol( Sxml.filedata.fname ); 
    seg8t.isCTseg = 0; 
    seg8t.help    = struct( ...
      'main'   , 'The main field "seg8t" includes single values from the SPM preprocessing used by SPM12, CAT12, or CTseg.', ...
      'Affine' , 'Affine transformation matrix from individual to (MNI) template space.', ...
      'lkp'    , 'Gaussians for each TPM class.', ...
      'wp'     , '', ...'proportion of each TPM class.', ...
      'mg'     , 'Weighting within each TPM class defined by lkp.', ...
      'mn'     , 'Gaussian peak value of each TPM class define by lkp.', ...
      'vr'     , 'Variance of each Gaussian peak value of each TPM class define by lkp.', ...
      'll'     , 'Final total log-likelyhood of the SPM preprocessing.', ...
      'isCTseg', 'Addditional flat that CTseg was used to process CT data.');
    vx_vol        = Sxml.qualitymeasures.res_vx_vol;
% ###################
% In case of CAT we could directly use the QC ratings.
% ###################
  end
  



  %% check for problems and skip in the worst case
  if max(seg8t.lkp) ~= 6 
    cat_io_cprintf('err','ERROR: Only 6 class models are supported!\n');
    return
  end
  % #####################
  % check number of Gaussian peaks per class ? 
  
  

 
  % == evaluate SPM mat entries to name some basic images properties ==
  % SPM main tissue thresholds
  tis.help.main       = ['The main "tis" structure inlcudes SPM-based measures (seg8*,WMth), ' ...
                         'image resolution (res_*) and major image class values (WM,GM,CSF,bone[cortex|marrow],head,background).'];
  tis.help.seg8o      = 'SPM seg8 main tissue intensity (mn of max mg in lkp)';
  tis.help.seg8ov     = 'SPM seg8 main tissue variance (vr of max mg in lkp)';
  tis.help.seg8n      = 'SPM main tissue intensity (mn of max mg in lkp) normalized by the WM'; % normalized by WM
  tis.help.seg8nv     = 'SPM main tissue variance (vr of max mg in lkp) normalized by the WM';
  tis.help.seg8con    = 'mimimum brain tissue contrast in SPM seg8t'; 
  tis.help.seg8conn   = 'mimimum brain tissue contrast in SPM seg8t normalized by the WM'; 
  tis.help.seg8CNR    = 'Noise estimate as minimum of the WM and CSF variance in percent (similar to BWP).';
  % ---
  tis.help.WMth       = 'SPM WM tissue intensity.'; 
  tis.help.res_vx_vol = 'Image voxel resolution in mm.'; 
  tis.help.res_RES    = 'RMS voxel resolution.'; 
  
  % create intensity variables
  tis.seg8o           = nan(1,6);
  tis.seg8ov          = nan(1,6);
  for ci = 1:max(seg8t.lkp) 
    % The SPM Gaussians seem to be unsorted and sorting based on the mean
    % value or the variance would be useful 
    sortvar = 'vr';
    [~,sorti] = sort( seg8t.(sortvar)(seg8t.lkp==ci) ); 
    var = seg8t.mn( seg8t.lkp==ci ); tis.seg8mns( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.mg( seg8t.lkp==ci ); tis.seg8mgs( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.vr( seg8t.lkp==ci ); tis.seg8vrs( seg8t.lkp==ci ) = var(sorti); 
    
    % How to average values ... well, when we are interested in changes of
    % intensities in bone(marrow) and the proportion of fat then (weighted)  
    % averaging should be fine (will depend on a subject and protocol).
    % To only use the strongest value is less stable (eg. for the 2-class CSF)
    tis.seg8o(ci)     = mean(seg8t.mg(seg8t.lkp==ci)' .* seg8t.mn(seg8t.lkp==ci),2);
    tis.seg8ov(ci)    = mean(seg8t.mg(seg8t.lkp==ci)' .* shiftdim(seg8t.vr(seg8t.lkp==ci),4),2);
  end

  % test if an image is a CT based on the intensities
  % low tissue contrast, but high bone, and low background
  if out.CTseg
    % test CT condition - use fixed values
    tis.WMth          = 40; 
    tis.seg8o         = [40 30 20 1024 0 -1024];
    tis.seg8n         = [40 30 20 1024 0 -1024] / 2000 +1;
    tis.seg8nv        = [40 30 20 1024 0 -1024];
    tis.seg8con       = nan; 
    tis.seg8conr      = nan; 
    tis.seg8CNR       = min( seg8t.vr( seg8t.lkp(:)==2 | seg8t.lkp(:)==3  | seg8t.lkp(:)==6) ) * 10e7; 
  else
    tis.WMth          = tis.seg8o(2);
    tis.seg8n         = tis.seg8o  ./ tis.WMth; % normalized by WM
    tis.seg8nv        = tis.seg8ov ./ tis.WMth;
    tis.seg8con       = min(abs(diff(tis.seg8o)));
    tis.seg8conr      = min(abs(diff(tis.seg8n)));
    tis.seg8CNR       = min( seg8t.vr(:) ./ tis.WMth ./ ...
      (seg8t.lkp(:)==2 | seg8t.lkp(:)==3  | seg8t.lkp(:)==6) ) .* tis.seg8conr * 3; 
  end

  % resolution measures
  tis.res_vx_vol    = vx_vol; 
  tis.res_RES       = mean(vx_vol.^2).^.5;
  
 
  % image weighting 
  tis.help.weighting = 'Image weighting based on SPM seg8t intensities (0=PDw; 1=T1w; 2=T2w; 3=MTw; 4=IRw, -1=CT).'; 
  if out.CTseg
      tis.weighting  = -1;
      tis.weightingn = 'CT';
  else
    if     tis.seg8n(3) < tis.seg8n(1)  &&  tis.seg8n(1) < tis.seg8n(2)      % T1: CSF < GM < WM 
      tis.weighting  = 1; 
      tis.weightingn = 'T1w';
    elseif tis.seg8n(3) > tis.seg8n(1)  &&  tis.seg8n(1) > tis.seg8n(2)      % T2: CSF > GM > WM
      tis.weighting  = 2;
      tis.weightingn = 'T2w';
    elseif tis.seg8o(3) < 0  && is.seg8o(2) < 2                              % MT: negative CSF values and not too high WM values
      tis.weighting  = 3;
      tis.weightingn = 'MTw';
    elseif tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) < .3                 % high int - low var
      tis.weighting  = 4;
      tis.weightingn = 'IRw'; % inverse recovery
    else
      tis.weighting  = 0; 
      tis.weightingn = 'PDw';
    end
  end


  % background types 
  tis.help.highBG = ['Intensity of the background based on the SPM seg8t intensities ' ...
    '(0=low,classical MRI; 1=high int low var, e.g. MP2Rage; 2=high int low var, e.g. IR; 3=mid int high var, e.g. MT)'];
  if     tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) > .3 % high int - high var
    tis.highBG     = 1;
    tis.highBGn    = 'high'; % intensity, high variance background (e.g. uncorrected MP2Rage)';
  elseif tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) < .3 % high int - low var
    tis.highBG     = 2;
    tis.highBGn    = 'high2'; % intensity, low variance background (e.g. inverse recovery)';
  elseif tis.seg8o(6) < tis.seg8o(3) 
    tis.highBG     = 0;
    tis.highBGn    = 'low'; % intensity, low variance background (e.g. classical MRI)';
  else
    tis.highBG     = 3;
    tis.highBGn    = 'mid'; % intensity, high variance background (e.g. MT)';
    tis.weighting  = 3;
    tis.weightingn = 'MTw';
  end


  % fat suppression
  tis.help.headFatType = 'Protocol intensity type of the head based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  tis.help.boneIntType = 'Protocol intensity type of the bone based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  minBone = min( seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) / tis.seg8o(2); % hard bone 
  maxBone = max( seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) / tis.seg8o(2); % bone marrow
  maxHead = max( seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) / tis.seg8o(2); % fat 
  if out.CTseg
    tis.headFatType   = 0; 
    tis.headFatTypen  = '-';
    tis.headBoneType  = 0; 
    tis.headBoneTypen = '-'; 
    % bone and head values of CTseg are not really useful :/ 
  else
    if maxHead > 1.2
      tis.headFatType  = 2; 
      tis.headFatTypen = 'high'; 
    elseif maxHead < tis.seg8n(3) 
      tis.headFatType  = 0; 
      tis.headFatTypen = 'low'; 
    else
      tis.headFatType  = 1; 
      tis.headFatTypen = 'mid'; 
    end
    if tis.headFatType > 1  &&  maxBone > 1
      tis.headBoneType  = 2; 
      tis.headBoneTypen = 'high'; 
    elseif tis.headFatType < 2  &&  minBone < 0.1
      tis.headBoneType  = 0; 
      tis.headBoneTypen = 'low'; 
    else
      tis.headBoneType  = 1; 
      tis.headBoneTypen = 'mid'; 
    end
  end  
   

  % skull-stripped ? 
  % #######################################################
  % create error / continue


  % defacing ?
  % deearing ?
  % >> affected area ? 
  %    surface percentage between cut (0-values) and uncut / full 
  %    background area
  % #######################################################
    
  
  % TIV estimate ? 
  % ############################
  % evaluate Affine value? > maybe with some fixed TPM TIV value?
  % evaluate log-likelihood?


  % final intensity measures
  tis.help.WM         = 'averaged SPM WM intensity';  % refined ( max in T1w, otherwise min )
  tis.help.GM         = 'averaged SPM GM intensity';  % no refinement possible
  tis.help.CSF        = 'averaged SPM CSF intensity'; % refined ( min in T1w, otherwise max )  
  tis.help.bonecortex = 'cortical bone intensity, i.e. min( seg8.mn( seg8t.lkp==4) )';  
  tis.help.bonemarrow = 'spongy bone intensity, i.e. max( seg8.mn( seg8t.lkp==4) )'; 
  tis.help.bonestruct = 'relation between cortical and spongy bone';
  tis.help.bone       = 'average bone intensity, i.e. mn( seg8.mn( seg8t.lkp==4 ) )';
  tis.help.head       = 'average head intensity, i.e. mn( seg8.mn( seg8t.lkp==5 ) )';
  tis.help.background = 'averaged SPM background intensity';
  % not so easy to extract the values from the typical 4 peaks 
  %tis.help.muscle     = 'normalized protocol-selected values of tis.*Head.'; 
  %tis.help.fat        = 'normalized protocol-selected values of tis.*Head.'; 
  
  tis.WM              = tis.seg8n(2);
  tis.GM              = tis.seg8n(1);
  tis.CSF             = tis.seg8n(3);
  tis.bonecortex      = min(seg8t.mn( seg8t.lkp==4 ) ) / tis.WMth;   
  tis.bonemarrow      = max(seg8t.mn( seg8t.lkp==4 ) ) / tis.WMth; 
  tis.bonedensity     = seg8t.mg( seg8t.mn == min(seg8t.mn( seg8t.lkp==4 ) ) )  + ... % percentage of min bone
                         0.5 * sum( seg8t.mg( seg8t.lkp==4 & ... % half percentage of median bone (neither min nor max) 
                          seg8t.mn ~= min(seg8t.mn( seg8t.lkp==4 ) ) & ...
                          seg8t.mn ~= max(seg8t.mn( seg8t.lkp==4 ) ) ));
  tis.bone            = tis.seg8n(4);
  tis.head            = tis.seg8n(5);        
  tis.background      = tis.seg8n(6);


  if verb == 3
    %% just an internal overview of SPM paraemter
    WMth = seg8t.mn( (seg8t.lkp==2) &  seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==2) )); % main WM intensity 
    fprintf('\nTissue table [class; propotion; WM normed intensity; WM normed var; normed var]:\n')
    disp( [seg8t.lkp; seg8t.mg';  
      seg8t.mn / WMth; shiftdim(seg8t.vr,1)  ./ WMth; % normalized by WM intensity 
      shiftdim(seg8t.vr,1) ./ WMth  .* (seg8t.mn / WMth);  ] ) 
     
  end

end
