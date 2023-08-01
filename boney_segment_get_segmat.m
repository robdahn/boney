function [ seg8t, tis, matm, mmatm, vx_vol ] = boney_segment_get_segmat(out,MAfn)   % TODO: (1) lkp tests; (2) skull-stripping/defacing tests; (3) TIV estimate  
%get_segmat. Get and evaluate SPM preprocessing structure.
%
%  [seg8t, tis, tis, matm, vx_vol] = get_spm8(P,MAfn)
%
%   P       .. one image
%   MAfn    .. fieldnames for matm
%
%   seg8t   .. spm8 structure without larger fields
%   tis     .. measures based on SPM tissue thresholds and 
%              basic information based on image and tissue properties
%   matm    .. line-wise output (print value)
%   mmatm   .. line-wise output (mark  value)
% 
% Comments: 
% * SPM tissue peaks for [GM WM CSF bone head background]
%   - default is [1 1 2 3 4 2]: 
%     . in T1 the lower CSF peak is the right one whereas the other one is PVE to GM or in the best case meninges?  
%       > would need a seperate class?
%     . the lower CSF threshold trend to overestimation in younger subjects - PVE?
%       > the CSF value is less robust for normalization 
%     . in mixed tissues, higher mg result often in higher vr that bias the mn  


% load SPM mat 
  seg8t          = load(out.P.seg8); 
  if out.CTseg
    seg8t.dat.model.gmm = rmfield(seg8t.dat.model.gmm,{'T','Sig'});
    % #### this is not fully working and the classe values are strange ...
    seg8t.lkp     = seg8t.sett.gmm.mg_ix;
    seg8t.mn      = seg8t.dat.model.gmm.m; 
    seg8t.mg      = seg8t.dat.model.gmm.gam'; 
    seg8t.vr      = seg8t.dat.model.gmm.W; 
     
    tmp           = spm_load_priors8(ps_fullfile(spm('dir'),'TPM','TPM.nii'));
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
  vx_vol         = sqrt(sum(seg8t.image(1).mat(1:3,1:3).^2));
  


  % check for problems and skip in worst case
  if numel(seg8t.tpm) ~= 6 
    cat_io_cprintf('err','ERROR: Only 6 class models are supported yet!\n');
    return
  end
  % #####################
  % check number of Gaussian peak per class ? 
  
  

 
  % == evaluate SPM mat entries to name some basic images properties ==
  % SPM main tissue thresholds
  tis.help.seg8o      = 'SPM seg8 main tissue intensity (mn of max mg in lkp)';
  tis.help.seg8ov     = 'SPM seg8 main tissue variance (vr of max mg in lkp)';
  tis.help.seg8n      = 'SPM main tissue intensity (mn of max mg in lkp) normalized by the WM'; % normalized by WM
  tis.help.seg8nv     = 'SPM main tissue variance (vr of max mg in lkp) normalized by the WM';
  tis.help.seg8con    = 'mimimum brain tissue contrast in SPM seg8t'; 
  tis.help.seg8conn   = 'mimimum brain tissue contrast in SPM seg8t normalized by the WM'; 
  tis.help.seg8CNR    = 'Noise estimate as minimum of the WM and CSF variance in percent (similar to BWP).';
  tis.help.WMth       = 'SPM WM tisse intensity.'; 
  tis.help.res_vx_vol = 'Image voxel resolution in mm.'; 
  tis.help.res_RES    = 'RMS voxel resolution.'; 
  tis.seg8o           = nan(1,6);
  tis.seg8ov          = nan(1,6);
  for ci = 1:max(seg8t.lkp) 
    % The SPM Gaussian's seams to be unsortet and sorting based on the mean
    % value or the variance would be useful 
    sortvar = 'vr';
    [~,sorti] = sort( seg8t.(sortvar)(seg8t.lkp==ci) ); 
    var = seg8t.mn( seg8t.lkp==ci ); tis.seg8mns( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.mg( seg8t.lkp==ci ); tis.seg8mgs( seg8t.lkp==ci ) = var(sorti); 
    var = seg8t.vr( seg8t.lkp==ci ); tis.seg8vrs( seg8t.lkp==ci ) = var(sorti); 
    
    % How to average values ... well, when we are interested in changes of
    % intensities in bone(marrow) and the propostion of fat then weighted  
    % averaging should be fine (will depend on subject and protocol).
    if 0 % use major class, ie. high proposion and low variance 
      tis.seg8o(ci)   = seg8t.mn( (seg8t.lkp==ci) & seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==ci)) );
      tis.seg8ov(ci)  = seg8t.vr( (seg8t.lkp==ci) & seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==ci)) );
    else % mix peak as average class intensity ... 
      tis.seg8o(ci)   = mean(seg8t.mg(seg8t.lkp==ci)' .* seg8t.mn(seg8t.lkp==ci),2);
      tis.seg8ov(ci)  = mean(seg8t.mg(seg8t.lkp==ci)' .* shiftdim(seg8t.vr(seg8t.lkp==ci),4),2);
    end
    if 0 %ci == 3 && sum(seg8t.lkp==ci)>1 % CSF subcase to avoid PVE in T1 images
      if mean( seg8t.mn(seg8t.lkp==ci) )  <  mean( seg8t.mn(seg8t.lkp<=2) ) % low CSF: CSF < WM+GM
        csfmn = seg8t.mn( seg8t.lkp==ci ); csfvr = seg8t.vr( seg8t.lkp==ci );
        [tis.seg8o(ci),csfmnmin] = min( csfmn );
        tis.seg8ov(ci) = csfvr(csfmnmin);
  
        tis.help.report_meninges = 'In case of multiple CSF Gaussians, the larger one in T1 often presents the intensity of meninges.'; 
        tis.meninges = mean( csfmn( setdiff(1:numel(csfmn),csfmnmin) ));
      else
        tis.help.report_meninges = 'In case of only one CSF Gaussian, the value of meninges should be in the middle between CSF and GM.'; 
        tis.meninges = mean( [tis.seg8ov(1) tis.seg8ov(3)] ); 
        % do not know enough about T2/PD case
      end
    end
  end

  % test if an image is a CT based on the intensities
  % low tissue contrast, but high bone, and low background
  %isCT = out.CTseg; % ###############  not working
       %  (tis.seg8o(1:3) > 0   & tis.seg8o(1:3) < 50) & ...  
       %  (tis.seg8o(4) > 600   & tis.seg8o(4) < 2000) & ...
       %  (tis.seg8o(5) > -50   & tis.seg8o(5) < 50)   & ...
       %  (tis.seg8o(6) > -2000 & tis.seg8o(6) < -600 );
 
  if out.CTseg
    % test CT condition 
    tis.WMth          = 40; 
    tis.seg8o         = [40 30 20 1024 0 -1024];
    tis.seg8n         = [40 30 20 1024 0 -1024] / 2000 +1;
    tis.seg8nv        = [40 30 20 1024 0 -1024];
    tis.seg8con       = nan; 
    tis.seg8conr      = nan; 
    tis.seg8CNR       = nan; 
  else
    tis.WMth          = tis.seg8o(2);
    tis.seg8n         = tis.seg8o  ./ tis.WMth; % normalized by WM
    tis.seg8nv        = tis.seg8ov ./ tis.WMth;
    tis.seg8con       = min(abs(diff(tis.seg8o)));
    tis.seg8conr      = min(abs(diff(tis.seg8n)));
    tis.seg8CNR       = min( seg8t.vr(:) ./ tis.WMth ./ (seg8t.lkp(:)==2 | seg8t.lkp(:)==3  | seg8t.lkp(:)==6) ) .* tis.seg8conr * 3; 
  end
  tis.res_vx_vol    = vx_vol; 
  tis.res_RES       = mean(vx_vol.^2).^.5;
  
  tis.spm_bone_con  = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
                      ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
            
  if max(seg8t.lkp) > 6
      tis.spm_bone_med1 = seg8t.mn( seg8t.lkp == 4 ) / tis.spm_bone_con;
      tis.spm_head_med1 = seg8t.mn( seg8t.lkp == 5 ) / tis.spm_bone_con;
  else
      % default case for 6 TPM classes
      tis.spm_con = sum(seg8t.mn( seg8t.lkp == 3 )' .* ...
        ( seg8t.mg( seg8t.lkp == 3 ) == max(seg8t.mg( seg8t.lkp == 3 ) ) ));
      [bone,bonei] = sort(seg8t.mn( seg8t.lkp == 4 ) / tis.spm_con); % order intensities within class 4 (normalized by CSF)
      bonevr = seg8t.vr( seg8t.lkp == 4 ) / tis.spm_con;
      if sum(seg8t.lkp == 4) == 3
        tis.spm_bone_med1   = bone(3); % highest value
        tis.spm_bone_med2   = bone(1); 
        tis.spm_bone_med3   = bone(2); 
        tis.spm_bone_medvr1 = bonevr(bonei(3)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(1)); 
        tis.spm_bone_medvr3 = bonevr(bonei(2)); 
      elseif sum(seg8t.lkp == 4) == 2 % CTseg
        tis.spm_bone_med1 = bone(1); % highest value
        tis.spm_bone_med2 = bone(2); 
        tis.spm_bone_medvr1 = bonevr(bonei(1)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(2)); 
      else
        tis.spm_bone_med1   = bone(1); % highest value
        tis.spm_bone_med2   = bone(1); 
        tis.spm_bone_medvr1 = bonevr(bonei(1)); % highest value
        tis.spm_bone_medvr2 = bonevr(bonei(1)); 
      end
  end


  % image weighting 
  tis.help.weighting = 'Image weigthing based on SPM seg8t intensities (0=PDw; 1=T1w; 2=T2w; 3=MTw; 4=IRw, -1=CT).'; 
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
    elseif tis.seg8o(3) < 0  && is.seg8o(2) < 2                              % MT: negative CSF values and not to high WM values
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
    '(0=low,classical MRI; 1=high int low var, eg. MP2Rage; 2=high int low var, eg. IR; 3=mid int high var, eg. MT)'];
  if     tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) > .3 % high int - high var
    tis.highBG     = 1;
    tis.highBGn    = 'high'; % intensity, high variance background (eg. uncorrected MP2Rage)';
  elseif tis.seg8o(6) > tis.seg8o(2) && tis.seg8nv(6) < .3 % high int - low var
    tis.highBG     = 2;
    tis.highBGn    = 'high2'; % intensity, low variance background (eg. inverse recovery)';
  elseif tis.seg8o(6) < tis.seg8o(3) 
    tis.highBG     = 0;
    tis.highBGn    = 'low'; % intensity, low variance background (eg. classical MRI)';
  else
    tis.highBG     = 3;
    tis.highBGn    = 'mid'; % intensity, high variance background (eg. MT)';
    tis.weighting  = 3;
    tis.weightingn = 'MTw';
  end


  % fat suppression
  tis.help.minBone     = 'Minimum Gaussian of the bone tissue class with more as 5% normalized by the WM intensity. ';
  tis.help.headFatType = 'Protocol intensity type of the head based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  tis.help.boneIntType = 'Protocol intensity type of the bone based on SPM seg8t values (0-low[<CSV], 1-mid[<WM], 2-[>WM]). ';
  tis.minHead = min(    seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.medHead = median( seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.maxHead = max(    seg8t.mn (seg8t.lkp==5 & seg8t.mg'>0.05)) ./ tis.seg8o(2); 
  tis.minBone = min(    seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  tis.medBone = median( seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  tis.maxBone = max(    seg8t.mn (seg8t.lkp==4 & seg8t.mg'>0.05)) ./ tis.seg8o(2);
  if out.CTseg
    tis.headFatType   = 0; 
    tis.headFatTypen  = '-';
    tis.boneIntType   = 0; 
    tis.boneIntTypen  = '-'; 
    % bone and head values of CTseg are not realy useful :/ 
  else
    if tis.maxHead > 1.2
      tis.headFatType  = 2; 
      tis.headFatTypen = 'high'; 
    elseif tis.maxHead < tis.seg8n(3) 
      tis.headFatType  = 0; 
      tis.headFatTypen = 'low'; 
    else
      tis.headFatType  = 1; 
      tis.headFatTypen = 'mid'; 
    end
    if tis.headFatType > 1  &&  tis.maxBone > 1
      tis.boneIntType   = 2; 
      tis.boneIntTypen  = 'high'; 
    elseif tis.headFatType < 2  &&  tis.minBone < 0.1
      tis.boneIntType   = 0; 
      tis.boneIntTypen  = 'low'; 
    else
      tis.boneIntType   = 1; 
      tis.boneIntTypen  = 'mid'; 
    end
  end  
    




  % skull-stripped ? 
  % #######################################################
  % create error / continue


  % defacing ?
  % deearing ?
  % >> affected area ? 
  %    surface percentage betwee cutted (0-values) and uncutted / full 
  %    background area
  % #######################################################
    
  
  % TIV estimate ? 
  % ############################
  % evaluate Affine value? > maybe with some fix TPM TIV value?
  % evaluate log-likelihood?



  if 0
  %% just an internal overview of SPM paraemter
    WMth = seg8t.mn( (seg8t.lkp==2) &  seg8t.mg' == max( seg8t.mg' .* (seg8t.lkp==2) )); % main WM intensity 
    fprintf('Tissue table [class; propotion; WM normed intensity; WM normed var; normed var]:\n')
    disp( [seg8t.lkp; seg8t.mg';  
      seg8t.mn / WMth; shiftdim(seg8t.vr,1)  ./ WMth; % normalized by WM intensity 
      shiftdim(seg8t.vr,1) ./ WMth  .* (seg8t.mn / WMth);  ] ) 
     
  end




  % report tissue types
  % - help
  tis.help.report        = 'Different measures used in the command line report.';
  tis.help.report_Tw     = 'tis.weightingn';
  tis.help.report_Thbg   = 'tis.highBGn'; 
  tis.help.report_Tfat   = 'tis.headFatTypen'; 
  tis.help.report_Tbone  = 'tis.boneFatTypen'; 
  tis.help.report_Tres   = 'tis.res_RES'; 
  tis.help.report_Tcnr   = 'tis.seg8CNR';
  tis.help.report_Pll    = 'seg8t.ll'; 
  tis.help.report_BG     = 'tis.seg8n(6)';
  tis.help.report_CSF    = 'tis.seg8n(3)';
  tis.help.report_GM     = 'tis.seg8n(1)';
  tis.help.report_WM     = 'tis.seg8n(2)';
  tis.help.report_muscle = 'normalized protocol-selected values of tis.*Head.'; 
  tis.help.report_fat    = 'normalized protocol-selected values of tis.*Head.'; 
  tis.help.report_bone   = 'tis.minBone';  
  tis.help.report_marrow = 'tis.maxBone'; 
  tis.help.report_MBI    = 'bone mean intensity (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2'; 
  
  % - measures
  tis.report.Tw     = tis.weightingn; 
  tis.report.Thbg   = tis.highBGn; 
  tis.report.Tskull = tis.seg8n(4); 
  tis.report.Thead  = tis.seg8n(5); 
  tis.report.Tfat   = tis.headFatTypen; 
  tis.report.Tbone  = tis.boneIntTypen; 
  tis.report.Tres   = tis.res_RES; 
  tis.report.Tcnr   = tis.seg8CNR;
  %tis.report.Pll    = seg8t.ll; 
  % brain tissue peaks 
  tis.report.BG     = tis.seg8n(6);
  tis.report.CSF    = tis.seg8n(3);
  tis.report.GM     = tis.seg8n(1);
  tis.report.WM     = tis.seg8n(2);
  % head tissues peaks 
  if tis.weighting == 1 % T1 
    musth = (seg8t.lkp==5) & (seg8t.mg'>0.05) & (seg8t.mn > tis.seg8o(3)) & (seg8t.mn < tis.seg8o(2)); % head mn value between CSF and WM
    fatth = (seg8t.lkp==5) & (seg8t.mg'>0.05) & (seg8t.mn > tis.seg8o(2)); 
    if sum(musth)>0,  tis.report.muscle = min(seg8t.mn(musth)) / tis.seg8o(2); else; tis.report.muscle = tis.medHead; end
    if sum(fatth)>0,  tis.report.fat    = max(seg8t.mn(fatth)) / tis.seg8o(2); else; tis.report.fat    = tis.maxHead; end
  else
     tis.report.muscle = tis.medHead; 
     tis.report.fat    = tis.minHead; 
  end
  % bone tissues peaks
  tis.report.bone   = tis.minBone;  
  tis.report.marrow = tis.maxBone; 
  % ###################### here we need to account the percentage values ;-)  
  if tis.boneIntType == 1
    tis.report.MBI    = (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2; 
  else
    tis.report.MBI    = 1 - ( (tis.maxBone - tis.report.bone) / (tis.report.fat - tis.report.bone) * 2); 
  end
  tis.report.MED      = nan;  
  tis.report.BT       = nan; 
  tis.report.HDT      = nan; 
  
  % create one line for the output table depending on the MAfn defintion
  for fni = 1:numel(MAfn)
    matm{1,fni}  = tis.report.(MAfn{fni});
    if isnumeric(tis.report.(MAfn{fni}))
      mmatm(1,fni) = tis.report.(MAfn{fni});
    else
      mmatm(1,fni) = 1; 
    end
  end

end
