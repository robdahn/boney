function [Vo, Yo, Ym, Yc, Yc4, Ybonemarrow, tismri, Si, St, Sti, Affine] = ... 
  boney_segment_fst(P,job,seg8t,tis,out,i,vx_vol) 
  
  Affine = seg8t.Affine; 
  Ym = []; Yc = {}; %Ybonemarrow = []; 
  tismri = struct();

  %% get bias corrected original image
  if 1
    Vo = spm_vol(P{i});
    Yo = single(spm_read_vols(Vo));

    Pc4  = out(i).P.cls{4}; %fullfile(out(i).P.orgpp,sprintf('c%d%s%s',4,out(i).P.orgff,out(i).P.ee));
    Vc4  = spm_vol(Pc4); 
    Yc4  = single(spm_read_vols(Vc4));
    [Ytiv,redR] = cat_vol_resize(Yc4,'reduceV',vx_vol,4,8,'max');
    Ytiv = cat_vol_morph(Ytiv>.5,'ldc',4) & Ytiv<.5; 
    Ytiv = cat_vol_morph(Ytiv>.5,'o',1); 
    Ytiv = cat_vol_resize(Ytiv,'dereduceV',redR);

    Ybonemarrow = single(Yo/tis.seg8o(3)) .* (Yc4>.5) / tis.seg8o(3);
    
    % get measures
    tismri.TIV     = cat_stat_nansum(Ytiv(:)>0.5) .* prod(vx_vol) / 1000;
    tismri.den     = [nan nan nan cat_stat_nansum(Yc4(:))     .* prod(vx_vol) / 1000 nan nan]; 
    tismri.vol     = [nan nan nan cat_stat_nansum(Yc4(:)>0.5) .* prod(vx_vol) / 1000 nan nan]; 
    tismri.volr    = tismri.vol ./ tismri.TIV;
    tismri.Tth     = [nan nan nan cat_stat_nansum(Yo(Yc4(:)>0.5)) / tis.seg8o(2) nan nan]; 
    tismri.Tsd     = [nan nan nan cat_stat_nanstd(Yo(Yc4(:)>0.5)) / tis.seg8o(2) nan nan]; 
    mn             = cat_stat_kmeans(Ybonemarrow(Yc4(:)>.5),3); % maximum value 
    tismri.iBone   = mn(2); %cat_stat_nanmean(Si.facevertexcdata);
    tismri.iBonemn = mn;
    [tismri.volmn3,tismri.volsd3,tismri.volvx3]  = cat_stat_kmeans( Ybonemarrow(Yc4(:)>.5) ,3); % median thickness to avoid low and high outliers
    [tismri.volmn2,tismri.volsd2,tismri.volvx2]  = cat_stat_kmeans( Ybonemarrow(Yc4(:)>.5) ,2); % median thickness to avoid low and high outliers
    tismri.volmd  = cat_stat_nanmedian( Ybonemarrow(Yc4(:)>.5));
    tismri.volmn  = cat_stat_nanmean( Ybonemarrow(Yc4(:)>.5) );
    tismri.volsd  = cat_stat_nanstd( Ybonemarrow(Yc4(:)>.5) );
    tismri.voliqr = iqr( Ybonemarrow(Yc4(:)>.5) );
%          nout(si).tis.seg8conr; 
   

    %% write output maps
    if job.output.writevol
      %%
      tdim = seg8t.tpm(1).dim; 
      M0   = seg8t.image.mat;          
      M1   = seg8t.tpm(1).mat;
  
      % affine and rigid parameters for registration 
      % if the rigid output is incorrect but affine is good than the Yy caused the problem (and probably another call of this function) 
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
      % midline map also for masking masking
      cat_io_writenii(Vo,Ybonemarrow,'',sprintf('bonemarrow%d_',job.opts.bmethod), ...
        'bone percentage position map','uint16',[0,0.001], ... 
        min([1 0 2],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.dartel]),trans);
      clear trans
    end
  end

  % no surfaces
  St = ''; Si = ''; Sti = '';