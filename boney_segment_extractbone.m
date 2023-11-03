function [Ybonepp, Ybonethick, Ybonemarrow, Yheadthick, vROI] = ...
  boney_segment_extractbone(Vo,Ym,Yc,Ye,Ya,Ymsk,seg8t,tis,out,job,vx_vol,YaROIname,RES,BB)
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

  %%
  if tis.weighting == -1
  %% CT images
    Ybrain0      = cat_vol_morph(cat_vol_morph((Yc{1} + Yc{2} + Yc{3})>.5,'lc',1,vx_vol),'lo',3,vx_vol); % remove 
    Ybraindist1  = cat_vbdist( single(Ybrain0>0.5) , (Yc{4} + Yc{5})>0 , vx_vol);
    Yhead        = Yc{1} + Yc{2} + Yc{3} + Yc{4}; Yhead(cat_vol_morph(Yhead>.5,'ldc',2)) = 1; 
    Ybone        = Yc{4};
    Ybrain       = (Yhead - Ybone) .* cat_vol_morph(Yhead>.5,'e') .* (Ybraindist1<4); % .* Ybrain;
    Ybraindist   = cat_vbdist( single(Ybrain>0.5) , Ybone>.5, vx_vol) .* (Ybone>0.5);
    Yheaddist    = cat_vbdist( single(Yhead<0.5)  , Ybone>.5, vx_vol) .* (Ybone>0.5);
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
    Ybraindist   = cat_vbdist( single(Ybrain>0.5) , Ybone1>0, vx_vol);
    Yheaddist    = cat_vbdist( single(Yhead>0.5)  , Ybone1>0, vx_vol);
    Ybonethick   = (Ybraindist + Yheaddist) .* (Ybone1>0);  % correct for voxel-size
    Ybonethick(Ybonethick > 100) = 0; 
    Ybonepp      = min(1,Yheaddist  ./ max(eps,Ybonethick)); Ybonepp(Ybrain>.5) = 1; % percentage map to
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
    
    % bonemarrow
    Ybonemarrow = single(Ym/tis.seg8n(3)) .* (Ybone>.5); % bias intensity normalized corrected
  end 


  % head 
  Yskull       = single(Ym/tis.seg8n(3)) .* (Yc{5}>.5);
  Ybgdist      = cat_vbdist( single(Yc{6}) , Ybrain<0.5, vx_vol);
  Ybndist      = cat_vbdist( cat_vol_morph(single(Ybrain + Ybone),'lc') , Yc{6}<.5, vx_vol);
  Yheadthick   = Ybndist + Ybgdist - max(0,Yheaddist - .5);
  Yheadthick   = cat_vol_localstat(Yheadthick,Yheadthick<1000,1,1);
  [~,YD]       = cat_vbdist(single(Yheadthick>0)); Yheadthick = Yheadthick(YD); % simple extension to avoid NaNs (BB limit is ok)
  clear Ybgdist Ybndist


 

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
  if ~isempty(job.opts.Pmask{1}), mskd = ' excluding the lower parts of the skull (masked)'; else, mskd = ' (unmasked)'; end 
  vROI.help = sprintf('ROI=0 is defined global values%s, whereas all other ROIs are without masking',mskd);
  for ri = 0:max(Ya(Ya(:)<intmax('uint16')))
    if ri == 0 || isnan(ri)
      % global values
      ri = 0; %#ok<FXSET> % case of failed atlas mapping given by NAN
      vROI.boneatlas_id(1,rii)      = inf;
      if ~isempty( job.opts.Pmask{1} ), vROI.boneatlas_name{1,rii} = 'full-masked'; 
      else,                             vROI.boneatlas_name{1,rii} = 'full-unmasked'; 
      end
      vROI.nonnanvol(1,rii)         = sum(Ya(:)>intmax('uint16')) ./ numel(Ya(:));
      vROI.bonemarrow(1,rii)        = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)>1  & Ybonemarrow(:)~=0 & Ybonepp(:)>.8 ) ); 
      vROI.bonecortex(1,rii)        = cat_stat_nanmean(   Ybonemarrow( Ymsk(:)>1  & Ybonemarrow(:)~=0 ) ); 
      vROI.bonethickness(1,rii)     = cat_stat_nanmean(   Ybonethick(  Ymsk(:)>1  & Ybonethick(:)~=0  ) ); 
      vROI.head(1,rii)              = cat_stat_nanmean(   Yskull(      Ymsk(:)>1  & Yskull(:)~=0      ) ); 
      vROI.headthickness(1,rii)     = cat_stat_nanmean(   Yheadthick(  Ymsk(:)>1  & Yheadthick(:)~=0  ) ); 
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
  [Ybonepp, Ybonethick, Ybonemarrow,Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceV' ,RES); % ############### INTERPOLATION ???
  [Ybonepp, Ybonethick, Ybonemarrow,Yheadthick] = cat_vol_resize({Ybonepp, Ybonethick, Ybonemarrow,Yheadthick} ,'dereduceBrain',BB); 

  if tis.headBoneType == 0 && tis.weighting > 0  &&  0 % #########
    % #### RD20231102: This correction is completely arbitrary and needs
    % explaination and further test (not working in OASIS3-1260 test-retest
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
    
    %%
    tdim = seg8t.tpm(1).dim; 
    M0   = seg8t.image.mat;          
    M1   = seg8t.tpm(1).mat;

    % affine and rigid parameters for registration 
    % if the rigid output is incorrect but affine is good then the Yy caused the problem (and probably another call of this function) 
    R               = spm_imatrix(seg8t.Affine); R(7:9)=1; R(10:12)=0; R=spm_matrix(R);  
    Mrigid          = M0\inv(R)*M1;                                                          % transformation from subject to registration space (rigid)
    Maffine         = M0\inv(seg8t.Affine)*M1;                                                 % individual to registration space (affine)
    mat0a           = seg8t.Affine\M1;                                                         % mat0 for affine output
    mat0r           = R\M1;                                                                  % mat0 for rigid ouput
    
    % settings for new boundary box of the output images 
    trans.native.Vo = seg8t.image(1); 
    trans.native.Vi = seg8t.image(1);
    trans.affine    = struct('odim',tdim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',seg8t.Affine);  % structure for cat_io_writenii
    trans.rigid     = struct('odim',tdim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);           % structure for cat_io_writenii
  
    job.output.bonemarrow  = struct('native',1,'warped',0,'dartel',3);
    job.output.position    = struct('native',0,'warped',0,'dartel',0);
    job.output.bone        = struct('native',0,'warped',0,'dartel',0);
    job.output.bonethick   = struct('native',0,'warped',0,'dartel',0);
    job.output.headthick   = struct('native',0,'warped',0,'dartel',0);

    % midline map also for masking masking
    %cat_io_writenii(Vo,Ybonemid,out.P.mridir,sprintf('bonemid%d_',job.opts.bmethod), ...
    %  'bone percentage position midline map','uint8',[0,1/255], ... 
    %  min([1 0 2],[job.output.position.native job.output.position.warped job.output.position.dartel]),trans); 
    % masked map for averaging
    Vo.fname  = fullfile(out.P.orgpp, [out.P.orgff, out.P.ee]);
    cat_io_writenii(Vo,Ybonepp,out.P.mridir,sprintf('bonepp%d_',job.opts.bmethod), ...
      'bone percentage position map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bone.native job.output.bone.warped job.output.bone.dartel]),trans);
    cat_io_writenii(Vo,Ybonemarrow,out.P.mridir,sprintf('bonemarrow%d_',job.opts.bmethod), ...
      'bone percentage position map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.dartel]),trans);
    cat_io_writenii(Vo,Ybonethick,out.P.mridir,sprintf('bonethickness%d_',job.opts.bmethod), ...
      'bone thickness map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bonethick.native job.output.bonethick.warped job.output.bonethick.dartel]),trans);
    cat_io_writenii(Vo,Yheadthick,out.P.mridir,sprintf('headthickness%d_',job.opts.bmethod), ...
      'head thickness map','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.headthick.native job.output.headthick.warped job.output.headthick.dartel]),trans);
  end
end
