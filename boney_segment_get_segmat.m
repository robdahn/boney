function [ seg8t, tis, vx_vol, trans] = boney_segment_get_segmat(out,job_ouput_writevol,verb)   
%boney_segment_get_segmat. Get and evaluate SPM preprocessing structure.
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
    seg8t           = load(out.P.seg8); 
    
    trans           = getDeformation(seg8t, job_ouput_writevol >= 3 ); 
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
    Sxml          = cat_io_xml(out.P.seg8); 
    if isfield(Sxml,'error')
      %tis    = struct();
      %vx_vol = [];
      cat_io_cprintf('err','boney.CATpperror: CAT prepreprocessing was not successful. Will went on with next subject. \n'); 
      %return
    end
    if isfield(Sxml,'SPMpreprocessing')
      seg8t       = Sxml.SPMpreprocessing;
    else
      seg8t       = struct('Affine',nan(4),'Affine0',nan(4),'lkp',nan(6),'mn',nan(6),'vr',nan(1,1,6),'ll',nan,'mg',nan(6,1));
    end
    tmp           = spm_load_priors8(fullfile(spm('dir'),'TPM','TPM.nii'));
    seg8t.tpmA    = tmp;
    seg8t.tpm     = rmfield(tmp.V,'private');
    if ~isfield(seg8t,'mg') % older cat versions don't have this field!
      for i = 1:max(seg8t.lkp)
        seg8t.mg(seg8t.lkp == i,1) = 1 / sum(seg8t.lkp == i); 
      end
    end
    if exist(Sxml.filedata.fname,'file')
      seg8t.image = spm_vol( Sxml.filedata.fname ); 
    else
      seg8t.image = spm_vol( Sxml.filedata.Fm ); 
      seg8t.image.fname = Sxml.filedata.fname; 
    end
    seg8t.isCTseg = 0; 
    seg8t.help    = struct( ...
      'main'   , 'The main field "seg8t" includes single values from the SPM preprocessing used by SPM12, CAT12, or CTseg.', ...
      'Affine' , 'Affine transformation matrix from individual to (MNI) template space.', ...
      'lkp'    , 'Gaussians for each TPM class.', ...
      'wp'     , '', ...'proportion of each TPM class.', ...
      'mg'     , 'Weighting within each TPM class defined by lkp.', ...
      'mn'     , 'Gaussian peak value of each TPM class define by lkp.', ...
      'vr'     , 'Variance of each Gaussian peak value of each TPM class define by lkp.', ...
      'll'     , 'Final total log-likelihood of the SPM preprocessing.', ...
      'isCTseg', 'Additional flat that CTseg was used to process CT data.');
    if isfield(Sxml,'qualitymeasures') && isfield(Sxml.qualitymeasures,'res_vx_vol') 
      vx_vol      = Sxml.qualitymeasures.res_vx_vol;
    else
      vx_vol      = nan(1,3); 
    end
% ###################
% In case of CAT we could directly use the QC ratings.
% ###################

    %% check for problems and skip in the worst case
    if ~isfield(Sxml,'error') && max(seg8t.lkp) ~= 6 
      tis = struct(); 
      cat_io_cprintf('err','ERROR:boney.SPMpperror: Only 6 class models are supported yet! Continue with next subject\n');
      return
    end

    % load Y? image
    if exist(out.P.y,'file')
      trans = struct();
    else
      cat_io_cprintf('err','ERROR:boney.CATpperror: No non-linear deformation available. Procced with affine data.\n');
      trans = struct(); 
    end
  end
  



  % #####################
  % check number of Gaussian peaks per class ? 
  
  

 
  % == evaluate SPM mat entries to name some basic images properties ==
  % SPM main tissue thresholds
  tis.help.main       = ['The main "tis" structure inlcudes SPM-based measures (seg8*,WMth), ' ...
                         'image resolution (res_*) and major image class values (WM,GM,CSF,bone[cortex|marrow],head,background).'];
  tis.help.seg8o      = 'SPM seg8 main tissue intensity (mn of max mg in lkp).';
  tis.help.seg8ov     = 'SPM seg8 main tissue variance (vr of max mg in lkp).';
  tis.help.seg8n      = 'SPM main tissue intensity (mn of max mg in lkp) normalized by the WM.'; % normalized by WM
  tis.help.seg8nv     = 'SPM main tissue variance (vr of max mg in lkp) normalized by the WM.';
  tis.help.seg8con    = 'Minimum brain tissue contrast in SPM seg8t.'; 
  tis.help.seg8conn   = 'Minimum brain tissue contrast in SPM seg8t normalized by the WM.'; 
  tis.help.seg8CNR    = 'Noise estimate as minimum of the WM and CSF variance in percent (similar to BWP).';
  % ---
  tis.help.WMth       = 'SPM WM tissue intensity.'; 
  tis.help.res_vx_vol = 'Image voxel resolution in mm.'; 
  tis.help.res_RES    = 'RMS voxel resolution.'; 
  
  % create intensity variables
  tis.seg8o           = nan(1,6);
  tis.seg8ov          = nan(1,6);
  if ~exist('Sxml','var') || ( exist('Sxml','var') && ~isfield(Sxml,'error') )
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
  tis.help.headFatType = 'Protocol intensity type of the head based on SPM seg8t values (0-low[<CSF], 1-mid[<WM], 2-[>WM]). ';
  tis.help.boneIntType = 'Protocol intensity type of the bone based on SPM seg8t values (0-low[<CSF], 1-mid[<WM], 2-[>WM]). ';
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
      % #### bug: not working in OASIS3-1260 (detect fat-suppression in on of the rescans
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
  tis.help.WM         = 'Averaged SPM WM intensity (default 1 class).';  % refined ( max in T1w, otherwise min )
  tis.help.GM         = 'Averaged SPM GM intensity (default 1 class).';  % no refinement possible
  tis.help.CSF        = 'Lower(T1)/higher(T2/PD) SPM CSF intensity (default 2 classes), average for more than 2 classes.'; 
                          % refined ( min in T1w, otherwise max )  
  tis.help.bonecortex = 'Cortical bone intensity, i.e. "min( seg8.mn( seg8t.lkp==4) )".';  
  tis.help.bonemarrow = 'Spongy bone intensity, i.e. "max( seg8.mn( seg8t.lkp==4) )".'; 
  tis.help.bone       = 'Weighted average of the SPM bone intensities, i.e. "mean(seg8t.mg(seg8t.lkp==4)'' .* seg8t.mn(seg8t.lkp==4),2)".';
  tis.help.head       = 'Weighted average of the SPM head intensities, i.e. "mean(seg8t.mg(seg8t.lkp==5)'' .* seg8t.mn(seg8t.lkp==5),2)".';
  tis.help.background = 'Averaged SPM background intensity.';
  % not so easy to extract the values from the typical 4 peaks 
  %tis.help.muscle     = 'Normalized protocol-selected values of tis.*Head.'; 
  %tis.help.fat        = 'Normalized protocol-selected values of tis.*Head.'; 
  
  tis.WM              = tis.seg8n(2);
  tis.GM              = tis.seg8n(1);
  if sum( seg8t.lkp == 3 ) == 2 &&  tis.weighting 
    %% take the value that stronger vary from the other tissues 
    tis.CSF = seg8t.mn( seg8t.lkp(:) == 3 & seg8t.mg(:) > .1) / tis.seg8o(2); 
    GMdiff  = abs(tis.CSF - tis.GM);
    tis.CSF = tis.CSF( GMdiff==max(GMdiff) ); 
  else
    tis.CSF           = tis.seg8n(3);
  end
  
  if 1 % no good correlation 
    [~,bminid]  = min( seg8t.mn  + inf*(seg8t.lkp==4) ); 
    if  tis.headFatType
      [~,hdmaxid] = max( seg8t.mn .* (seg8t.lkp==5) ); 
    else
      [~,hdmaxid] = min( seg8t.mn + inf*(seg8t.lkp==5 | seg8t.mn > min( seg8t.mn + inf*(seg8t.lkp==5) )  ) ); 
    end
    tis.bonecortexpercentage  = seg8t.mg( bminid ); 
    tis.headfatpercentage     = seg8t.mg( hdmaxid ); 
  end
  tis.bonecortex            = min(seg8t.mn( seg8t.lkp==4 ) ) / tis.WMth;   
  tis.bonemarrow            = max(seg8t.mn( seg8t.lkp==4 ) ) / tis.WMth; 
  tis.bone                  = tis.seg8n(4);
  tis.head                  = tis.seg8n(5);        
  tis.background            = tis.seg8n(6);

  % normalization values as minimal main CSF or background and maximal WM intensity
  tis.intnorm = [
    min([ seg8t.mn(seg8t.lkp(:) == 3              & seg8t.mg(:) > .1) , ...
          seg8t.mn(seg8t.lkp(:) == max(seg8t.lkp) & seg8t.mg(:) > .1) ] ), ...
    max(  seg8t.mn(seg8t.lkp(:) == 2              & seg8t.mg(:) > .3)   ) ];


  if verb == 3
    %% just an internal overview of SPM paraemter
    WMth = seg8t.mn( (seg8t.lkp==2) &  seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==2) )); % main WM intensity 
    fprintf('\nTissue table [class; proportion; WM normed intensity; WM normed var; normed var]:\n')
    disp( [seg8t.lkp; seg8t.mg';  
      seg8t.mn / WMth; shiftdim(seg8t.vr,1)  ./ WMth; % normalized by WM intensity 
      shiftdim(seg8t.vr,1) ./ WMth  .* (seg8t.mn / WMth);  ] ) 
     
  end

end
function trans = getDeformation(seg8t, dodefs)
%% relevant lines from spm_preproc_write8 and cat_main_registration

  tdim = seg8t.tpm(1).dim; 
  M0   = seg8t.image.mat;          
  M1   = seg8t.tpm(1).mat;

  % affine and rigid parameters for registration 
  % if the rigid output is incorrect but affine is good then the Yy caused the problem (and probably another call of this function) 
  R               = spm_imatrix(seg8t.Affine); R(7:9)=1; R(10:12)=0; R=spm_matrix(R);  
  Mrigid          = M0\inv(R)*M1;                      % transformation from subject to registration space (rigid)
  Maffine         = M0\inv(seg8t.Affine)*M1;           % individual to registration space (affine)
  mat0a           = seg8t.Affine\M1;                   % mat0 for affine output
  mat0r           = R\M1;                              % mat0 for rigid ouput
  M2              = inv(Maffine);                      % warped - if we use the US than we may have to use rigid
 
  % Use SPM deformation fir the actual dimensions and orientations of the tissue priors.
  if dodefs
    prm       = [3 3 3 0 0 0];
    Coef      = cell(1,3);
    Coef{1}   = spm_bsplinc(seg8t.Twarp(:,:,:,1),prm);
    Coef{2}   = spm_bsplinc(seg8t.Twarp(:,:,:,2),prm);
    Coef{3}   = spm_bsplinc(seg8t.Twarp(:,:,:,3),prm);
    d         = seg8t.image(1).dim(1:3);
    [x1,x2,~] = ndgrid(1:d(1),1:d(2),1); x3 = 1:d(3);
    M         = M1 \ seg8t.Affine * seg8t.image(1).mat;
    y         = zeros([seg8t.image(1).dim(1:3),3],'single'); % native < template
    for z=1:length(x3)
      [t1,t2,t3] = defs(Coef,z,seg8t.MT,prm,x1,x2,x3,M);
      y(:,:,z,1) = t1;
      y(:,:,z,2) = t2;
      y(:,:,z,3) = t3;
    end
  
    % inverse deformation
    if 0
      yi = spm_diffeo('invdef',y,tdim,eye(4),M0); % template > native (but in subject space)
      yi = spm_extrapolate_def(yi,M1);            % apply BB
      w  = max( eps , abs(spm_diffeo('def2det', yi ) ) ); 
      yi2 = yi; 
    else  
      yid             = spm_diffeo('invdef', y  , d, eye(4), M1\seg8t.Affine*M0); 
      yi              = spm_diffeo('invdef', yid, d, inv(M1\seg8t.Affine*M0), eye(4)); clear yid; 
      yi2             = spm_diffeo('invdef', yi , tdim, eye(4), eye(4)); 
      w               = max( eps , abs(spm_diffeo('def2det', yi2 ) ) ); 
    end
  end

  % final structure
  trans.native.Vo = seg8t.image(1); 
  trans.native.Vi = seg8t.image(1);
  trans.affine    = struct('odim',tdim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',seg8t.Affine);  % structure for cat_io_writenii
  trans.rigid     = struct('odim',tdim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);           % structure for cat_io_writenii
  if dodefs
    trans.warped  = struct('y',y ,'yi',yi2,'w',w,'odim',tdim,'M0',M0,'M1',M1,'M2',M2,'dartel',1,'fs',0);  % nicer version with interpolation
  end
end
%==========================================================================
% function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
%==========================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
  iMT = inv(MT);
  x1  = x0*iMT(1,1)+iMT(1,4);
  y1  = y0*iMT(2,2)+iMT(2,4);
  z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
  
  % Eliminate NaNs (see email from Pratik on 01/09/23)
  x1  = min(max(x1,1),size(sol{1},1));
  y1  = min(max(y1,1),size(sol{1},2));
  z1  = min(max(z1,1),size(sol{1},3));
  
  x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
  y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
  z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
  x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
  y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
  z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
end


