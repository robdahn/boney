function [Ybonepp, Ybonethick, Ybonemarrow, Yheadthick, vROI] = ...
  boney_segment_extractbone(Vo,Ym,Yc,Ye,Ya,Ymsk,trans,seg8t,tis,tismri,out,job,vx_vol,YaROIname,RES,BB)
%% * Report: 
%   - better an upper slice?
%   - optimize print font size
%   - add basic parameters (reduce,mask,Atlas,tpm)
%   - add tissue volumes
%   - add tissue intensities
%   - add histograms for tissues
%   - modify colorbar position and labeling > CAT
%   - use fat/muscle colors (yellow, pink)  
%   - use green/cyan for bone?
%   - affine registration surf problem
%   - use final segmentation for overlay but mark outliers 
%   - Opimize report line >> table like similar to QC with vols, thick & intensities 

  nvbdist      = 3; 

  %%
  if tis.weighting == -1
  %% CT images
    Ybrain0      = cat_vol_morph(cat_vol_morph((Yc{1} + Yc{2} + Yc{3})>.5,'lc',1,vx_vol),'lo',3,vx_vol); % remove 
    Ybraindist1  = cat_vbdist( single(Ybrain0>0.5) , (Yc{4} + Yc{5})>0 , vx_vol);
    Yhead        = Yc{1} + Yc{2} + Yc{3} + Yc{4}; Yhead(cat_vol_morph(Yhead>.5,'ldc',2)) = 1; 
    Ybone        = Yc{4};
    Ybrain       = (Yhead - Ybone) .* cat_vol_morph(Yhead>.5,'e') .* (Ybraindist1<4); % .* Ybrain;
    Ybraindist   = vbx_dist( single(Ybrain>0.5) , Ybone>.5, vx_vol, 1) .* (Ybone>0.5);
    Yheaddist    = vbx_dist( single(Yhead<0.5)  , Ybone>.5, vx_vol, 1) .* (Ybone>0.5);
    Ybonethick   = Ybraindist  + Yheaddist;  % correct for voxel-size
    Ybonepp      = min(1,Yheaddist  ./ max(eps,Ybonethick));  Ybonepp(Ybrain>.5) = 1; % percentage map to
    Ybonemarrow  = Ym;

  else
  % MRI images
    if 0
      Yheadbone  = cat_vol_morph(cat_vol_morph( smooth3(Yc{4} + Yc{5} + Yc{6}) > 0.5,'c',3),'e');     % create 
      Yheadbone  = cat_vol_smooth3X(Yheadbone,4)>.5;
      Ybrain     = single(Yc{1} + Yc{2} + Yc{3} - Yheadbone);                                 % SPM brain 
      Ybrainc    = (Yc{1} + Yc{2} + Yc{3}) .* (Ybrain<.5 & Yheadbone); 
      Yc{4}      = Yc{4} + Ybrainc; 
      Ybrain     = max( Ybrain , cat_vol_smooth3X(Ybrain,4)>.5); % remove vein
      %%
      Yhead      = min(1,single(Yc{5} + Yc{6}));
      Yhead(~cat_vol_morph(cat_vol_morph(Yhead>.5,'lo',1),'d') & Yhead>.5 & Ybrain<.5) = 0.5; 
      Ybg        = single(Yc{6});
      Ybone      = single(Yc{4});
     
      %% 
      if 1
        Ybrain = single(Ybrain); Ybone = single(Ybone); Yhead = single(Yhead);
        spm_smooth(Ybrain,Ybrain,2./vx_vol);
        spm_smooth(Yhead,Yhead,2./vx_vol);
        spm_smooth(Ybg,Ybg,2./vx_vol);
      end
      Ybone     = ~(Ybrain>.5 | Yhead>.5); 
    else
      Ybrain     = single(Yc{1} + Yc{2} + Yc{3});
      Yhead      = single(Yc{5} + Yc{6});
      Ybone      = single(Yc{4});
    end
    if ~isempty(Ye)
      Ybrain = Ybrain + Ye{1};
    end
    
    % remove PVE
    %  Ybonee     = Ybone>.5 & ~cat_vol_morph(Ybone>.5,'e');
    %  Ybonem     = cat_vol_localstat(single(Yo + ((Ybone<.5)*100)),Ybone>.5,2,2); 
    
    Ybrain = cat_vol_morph( cat_vol_morph( Ybrain>.5 , 'lc' , 2) , 'lo' , 2);
    Yhead  = ~cat_vol_morph( Yhead<.5 , 'lo' , 1); 

    Ybrain = cat_vol_smooth3X( Ybrain , 1)>.5; 
    Yhead  = cat_vol_smooth3X( Yhead  , 1)>.5; 
    Ybone1 = 1 - Ybrain - Yhead; 
    
  

    
    %% bone layers
    [Ybrainr,Yheadr,Ybone1r,res] = cat_vol_resize({Ybrain,Yhead,Ybone1},'reduceV',vx_vol,2,64,'meanm');
    Ybraindist   = vbx_dist( Ybrainr , Ybone1r>0, res.vx_volr, nvbdist, 0);
    Yheaddistr   = vbx_dist( Yheadr  , Ybone1r>0, res.vx_volr, nvbdist, 0);
    Ybonethick   = (Ybraindist + Yheaddistr) .* (Ybone1r>0);  % correct for voxel-size
    Ybonethick   = cat_vol_approx(Ybonethick,'rec'); 
    Ybonethick   = cat_vol_smooth3X(Ybonethick,1); 
    Ybonepp      = min(1,Yheaddistr ./ max(eps,Ybonethick)); Ybonepp(Ybrainr>.5) = 1; % percentage map to
    Ybonethick   = cat_vol_resize(Ybonethick,'dereduceV',res);
    Ybonepp      = cat_vol_resize(Ybonepp,'dereduceV',res);
    clear Ybrainr Yheadr Ybone1r Ybraindist; 
    %%
    if 0
      % head values
      Ybonehead = ~(Ybrain>.5 | Ybg>.5); 
      %Ybonehead    = Ybonehead .* (Ybraindist<30 & Ybonehead); 
      % Ybg          = max(Ybg,Ybraindist>=30 & Ybonehead); 
      Ybraindist2  = cat_vbdist( single(Ybrain>0.5) , Ybonehead, vx_vol) .* (Ybonehead>0.5);
      Yheaddist2   = cat_vbdist( single(Ybg>0.5)    , Ybonehead, vx_vol) .* (Ybonehead>0.5);
      Ybonethick2  = Ybraindist2 + Yheaddist2;  % correct for voxel-size
      Ybonepp2     = min(1,Yheaddist2 ./ max(eps,Ybonethick2)); Ybonepp2(Ybrain>.5) = 1;    % percentage map to
      Ybonepp      = max(Ybonepp,max(0,Ybonepp2*3-2));
    end
    clear braindist Ybonedist
    
  end 


  %% head 
% PK20240405: check normalization for skull - using the CSF value may lead to unwanted statistical dependencies ! 
  Yskull       = single(Ym/tis.seg8n(3)) .* (Yc{5}>.5); 
  [Yc6,Yheadr,Ybrainr,res] = cat_vol_resize({single(Yc{6}),cat_vol_morph(single(Ybrain + Ybone),'lc'),Ybrain},'reduceV',vx_vol,2,64,'meanm');
  Ybgdist      = vbx_dist( Yc6 , Ybrainr<0.5, res.vx_volr, nvbdist,0);
  Ybndist      = vbx_dist( Yheadr , Yc6<.5, res.vx_volr, nvbdist,0);
  Yheadthick   = (Ybndist + Ybgdist - Yheaddistr) .* (Ybrainr<0.5 & Yc6<.5); 
  clear Yheadr Yc6 Ybrainr Ybgdist Ybndist Yheaddistr; 
  % headthickness
  Yheadthick   = cat_vol_approx(Yheadthick,'rec');
  Yheadthick   = cat_vol_smooth3X(Yheadthick,1); 
  Yheadthick   = cat_vol_resize(Yheadthick,'dereduceV',res); clear res; 

%%

  % define normalized bonemarrow 
  %Ym = (Yo - min( tis.seg8o )) / (tis.seg8o(2) - min([tis.seg8o(3),tis.seg8o(end)]));
  Ybonemarrow = single( Ym ) .* Ybone;
  if tis.weighting >= 0
    % bonemarrow 
    switch job.bnorm
      case 'WM'
      case 'CSF' 
        Ybonemarrow = single( Ybonemarrow / tismri.CSF ); % bias intensity normalized corrected
      case 'muscle' 
        Ybonemarrow = single( Ybonemarrow / ( tismri.int.head_muscle / 4 ) ); % 
      case 'bone' 
        Ybonemarrow = single( Ybonemarrow / tismri.int.bone_cortex ); % bias intensity normalized corrected
      case 'fat'
        Ybonemarrow = single( Ybonemarrow / ( tismri.int.head_fat / 16 ) ); % 
    end
  end
 

%% ###################
% edge-measure ! 
% * this one is not realy working (in low fat cases?)
% * bone / bone-marrow segmentation: 
%   - bone marrow as outstanding maximum structure in the middle/center area
%   


  %% measures as column elements
  %  - first tested showed that mean/median works best for BMD and that SD/IQR are much worse 
  %  - we finally decided to keep only the mean as (i) it is more expected 
  %    and (ii) we already did some outlier correction by masking  
  rii = 1;
  vROI.help = [
      'A masked image is used (if available) for global values to extract only the upper part of the skull, ' ...
      'whereas no masking is used in case of atlas regions. '];
  for ri = 0:max(Ya(Ya(:)<intmax('uint16')))
    if ri == 0 || isnan(ri)
      % global values
      ri = 0; %#ok<FXSET> % case of failed atlas mapping given by NAN
      vROI.boneatlas_id(1,rii)      = inf;
      if ~isempty( job.opts.Pmask{1} ), vROI.boneatlas_name{1,rii} = 'full-masked'; 
      else,                             vROI.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      vROI.nonnanvol(1,rii)         = sum(Ya(:)>intmax('uint16')) ./ numel(Ya(:));
      vROI.bonemarrow(1,rii)        = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)  & Ybonemarrow(:)~=0 & Ybonepp(:)>.8 ) ); 
      vROI.bonecortex(1,rii)        = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)  & Ybonemarrow(:)~=0 ) ); 
      vROI.bonethickness(1,rii)     = cat_stat_nanmean(   Ybonethick(  Ymsk(:)  & Ybonethick(:)~=0  ) ); 
      vROI.head(1,rii)              = cat_stat_nanmean(   Yskull(      Ymsk(:)  & Yskull(:)~=0      ) ); 
      vROI.headthickness(1,rii)     = cat_stat_nanmean(   Yheadthick(  Ymsk(:)  & Yheadthick(:)~=0  ) ); 
      rii = rii + 1;
    else
      % regional values
      if sum(Ya(:)==ri)~=0
        vROI.boneatlas_id(1,rii)    = ri;  
        if isempty(YaROIname) %|| numel(YaROIname)>max(Ya(Ya(:)<intmax('uint16')))
          vROI.boneatlas_name{1,rii}  = sprintf('ROI%d',ri); 
        else
          vROI.boneatlas_name{1,rii}  = YaROIname{rii};
        end
        vROI.nonnanvol(1,rii)       = sum(Ya(:)==ri) ./ numel(Ya(:));
        vROI.bonemarrow(1,rii)      = cat_stat_nanmean(  Ybonemarrow( Ybonemarrow(:)~=0 & Ya(:)==ri & Ybonepp(:)>.8 ) ); 
        vROI.bonecortex(1,rii)      = cat_stat_nanmean(  Ybonemarrow( Ybonemarrow(:)~=0 & Ya(:)==ri) ); 
        vROI.bonethickness(1,rii)   = cat_stat_nanmean(  Ybonethick(  Ybonethick(:)~=0  & Ya(:)==ri) );
        vROI.head(1,rii)            = cat_stat_nanmean(  Yskull(      Yskull(:)~=0      & Ya(:)==ri) ); 
        vROI.headthickness(1,rii)   = cat_stat_nanmean(  Yheadthick(  Yheadthick(:)~=0  & Ya(:)==ri) ); 
        rii = rii + 1;
      end
    end
  end



  %% restore resolution & boundary box
  [Ybonepp, Ybonethick, Ybonemarrow, Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceV' ,RES); % ############### INTERPOLATION ???
  [Ybonepp, Ybonethick, Ybonemarrow, Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceBrain',BB); 

  if tis.headBoneType == 0 && tis.weighting > 0  &&  0 % #########
    % #### RD20231102: This correction is completely arbitrary and needs
    % explanation and further test (not working in OASIS3-1260 test-retest
    Ybonemarrow = Ybonemarrow * 3;
  end
  if tis.weighting == -1
    Ybonemarrow(Ybonemarrow==0) = -1024;
  end

  


  %% write output maps
  %  - what maps do we really need?
  %  - intensity normalized maps used for normalization and analysis? 
  %  - 
  if job.output.writevol
    
    %% see also boney_segment_get_seg8mat
    if ~exist('trans','var')
      tdim = seg8t.tpm(1).dim; 
      M0   = seg8t.image.mat;          
      M1   = seg8t.tpm(1).mat;
  
      % affine and rigid parameters for registration 
      % if the rigid output is incorrect, but the affine is good, then the Yy caused the problem (and probably another call of this function) 
      R               = spm_imatrix(seg8t.Affine); R(7:9)=1; R(10:12)=0; R=spm_matrix(R);  
      Mrigid          = M0\inv(R)*M1;                                                            % transformation from subject to registration space (rigid)
      Maffine         = M0\inv(seg8t.Affine)*M1;                                                 % individual to registration space (affine)
      mat0a           = seg8t.Affine\M1;                                                         % mat0 for affine output
      mat0r           = R\M1;                                                                    % mat0 for rigid ouput
      
      % settings for new boundary box of the output images 
      trans.native.Vo = seg8t.image(1); 
      trans.native.Vi = seg8t.image(1);
      trans.affine    = struct('odim',tdim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',seg8t.Affine);  % structure for cat_io_writenii
      trans.rigid     = struct('odim',tdim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);             % structure for cat_io_writenii
    end

    %% manual defintion of what volumes we write
    native = job.output.writevol==1 || job.output.writevol==4;
    warped = job.output.writevol==2 || job.output.writevol==4; 
    dartel = job.output.writevol==3 || job.output.writevol==4; 
    job.output.bonemarrow  = struct('native',native,'warped',warped,'mod',warped,'dartel',2*dartel); % dartel * 3 for affine & rigid
    job.output.position    = struct('native',0,'warped',0,'dartel',0);
    job.output.bonethick   = struct('native',0,'warped',0,'dartel',0);
    job.output.headthick   = struct('native',0,'warped',0,'dartel',0);

    % midline map also for masking masking
    %cat_io_writenii(Vo,Ybonemid,out.P.mridir,sprintf('bonemid%d_',job.opts.bmethod), ...
    %  'bone percentage position midline map','uint8',[0,1/255], ... 
    %  min([1 0 2],[job.output.position.native job.output.position.warped job.output.position.dartel]),trans); 
    % masked map for averaging
    Vo.fname  = fullfile(out.P.orgpp, [out.P.orgff, out.P.ee]);
    cat_io_writenii(Vo,Ybonemarrow,out.P.mridir,sprintf('bonemarrow%d_',job.opts.bmethod), ...
      'bone(marrow) map','uint16',[0,0.001], ... 
      min([1 1 1 3],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.mod job.output.bonemarrow.dartel]),trans);
    % no modulation here 
    cat_io_writenii(Vo,Ybonepp,out.P.mridir,sprintf('bonepp%d_',job.opts.bmethod), ...
      'bone percentage position map','uint16',[0,0.001], ... 
      min([1 1 3],[job.output.position.native job.output.position.warped job.output.position.dartel]),trans);
    cat_io_writenii(Vo,Ybonethick,out.P.mridir,sprintf('bonethickness%d_',job.opts.bmethod), ...
      'bone thickness map','uint16',[0,0.001], ... 
      min([1 1 3],[job.output.bonethick.native job.output.bonethick.warped job.output.bonethick.dartel]),trans);
    cat_io_writenii(Vo,Yheadthick,out.P.mridir,sprintf('headthickness%d_',job.opts.bmethod), ...
      'head thickness map','uint16',[0,0.001], ... 
      min([1 1 3],[job.output.headthick.native job.output.headthick.warped job.output.headthick.dartel]),trans);
  end
end
function Yd = vbx_dist(Yb,Ym,vx_vol,steps,sm)
  if ~exist('vx_vol','var'); vx_vol = [1 1 1]; end
  if ~exist('steps','var'); steps = 3; else, steps = round(steps); end
  if ~exist('sm','var'); sm = 0; end

  if sm, Yb = cat_vol_smooth3X(Yb,sm ./ vx_vol); end
  
  Yd = zeros(size(Yb),'single');
  for ss = 1/(steps+1):1/(steps+1):(1-1/(steps+1))
    Yd = Yd + 1/steps .* max(0,cat_vbdist( single( Yb > ss ) , Ym , vx_vol) - 0.5);
  end

  maxd = max( size(Ym)./vx_vol); 
  Yd(Yd>maxd | isnan(Yd) | isinf(Yd)) = 0;
  Yd(Yd<0) = 0; 
end