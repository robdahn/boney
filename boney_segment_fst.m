function [Vo, Ym, Yc, Ybonemarrow, tismri, Si, Stm, Affine] = ... 
  boney_segment_fst(Pi, Pc4, job, seg8t, tis, vx_vol) 
%boney_segment_fst. Fast bone processing, export and report preperation.  
% 
%  [Vo, Yo, Ym, Yc, Yc4, Ybonemarrow, tismri, Si, St, Sti, Affine] = ... 
%    boney_segment_fst(Pi, Pc4, job, seg8t, tis, vx_vol) 
%
%  Pi           .. ith input file
%  Pc4          .. bone class image
%  job          .. main job structure to include some parameters in the
%                  filename etc.
%  seg8t        .. spm8 structure without larger fields
%  tis          .. measures based on SPM tissue thresholds and basic 
%                  information based on image and tissue properties
%  vx_vol       .. voxel-size
%  Vo           .. original image header
%  Ym           .. intensity normalized image
%  Yc           .. segment class images (cell)
%  Ybonemarrow  .. bone marrow map (masked normalized image)
%  tismri       .. structure with MRI based measures
%  Si           .. bone intensity surface
%  Stm          .. bone thickness surface with regional boundaries%  
%  Affine       .. (reprocessed) affine transformation to MNI
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


% TODO: 
%  * maybe extract fast head measure?

  Affine = seg8t.Affine;      % this was reprocessed in other cases
  tismri = struct();          % MRI-base result structure 
  Yc = {};                    % no general tissue classes 
  Si = ''; Stm = '';          % no surfaces at all

  Vo = spm_vol(Pi);
  Yo = single(spm_read_vols(Vo));

  Vc4   = spm_vol(Pc4); 
  Yc{4} = single(spm_read_vols(Vc4));

  % approximation of TIV as enclosed volume
  [Ytiv,redR] = cat_vol_resize(Yc{4},'reduceV',vx_vol,4,8,'max');
  Ytiv = cat_vol_morph(Ytiv>.5,'ldc',4) & Ytiv<.5; 
  Ytiv = cat_vol_morph(Ytiv>.5,'o',1); 
  Ytiv = cat_vol_resize(Ytiv,'dereduceV',redR);

  % the bone intensity map
  Ybonemarrow          = single(Yo/tis.seg8o(3)) .* (Yc{4}>.5);
  Ym                   = Yo / tis.seg8o(2);
  
  % prepare measures just for the bone (class 4)
  tismri.help.TIV      = 'Total intracranial volume (GM + WM + CSF) - fast approximation.';
  tismri.help.vol      = 'Volume of the SPM tissues classes in mm (probability >.5).';
  tismri.help.volr     = 'Relative volume of the SPM tissues classes (probability >.5) normalized by TIV.';
  tismri.help.Tth      = 'Median intensity of the (optimized) tissue class (~peak intensity).'; 
  tismri.help.Tsd      = 'SD of the intensity of the (optimized) tissue class (~peak intensity).'; 
  tismri.help.bone_med = 'Median intensity of the SPM bone segment.'; 
  tismri.help.int      = ['Intensity based evaluated tissue classes: ' ...
    '(1) GM:   nan; ' ...
    '(2) WM:   nan; ' ...
    '(3) CSF:  nan; ' ...
    '(4) bone: kmeans with 3 classes for bone (1) and bone marrow (2); ' ...
    '(5) head: nan; ' ...
    '(6) bg:   nan. ' ...
    ];

  tismri.TIV          = cat_stat_nansum(Ytiv(:)>0.5) .* prod(vx_vol) / 1000;
  tismri.vol          = [nan nan nan cat_stat_nansum(Yc{4}(:)>0.5) .* prod(vx_vol) / 1000 nan nan]; 
  tismri.volr         = tismri.vol ./ tismri.TIV;
  tismri.Tth          = [nan nan nan cat_stat_nansum(Yo(Yc{4}(:)>0.5)) / tis.seg8o(2) nan nan]; 
  tismri.Tsd          = [nan nan nan cat_stat_nanstd(Yo(Yc{4}(:)>0.5)) / tis.seg8o(2) nan nan]; 
  tismri.bone_med     = cat_stat_nanmedian( Ybonemarrow(Yc{4}(:)>.5) );
  tismri.int          = struct('GM', nan, 'WM', nan, 'CSF', nan,...
                         'bone', cat_stat_kmeans(Ybonemarrow(Yc{4}(:)>.5),2), 'head', nan, 'BG', nan);


  %% write output maps
  if job.output.writevol
    tdim = seg8t.tpm(1).dim; 
    M0   = seg8t.image.mat;          
    M1   = seg8t.tpm(1).mat;

    % affine and rigid parameters for registration 
    % if the rigid output is incorrect but affine is good than the Yy  
    % caused the problem (and probably another call of this function) 
    R               = spm_imatrix(seg8t.Affine); R(7:9)=1; R(10:12)=0; R=spm_matrix(R);  
    Mrigid          = M0\inv(R)*M1;               % transformation from subject to registration space (rigid)
    Maffine         = M0\inv(seg8t.Affine)*M1;    % individual to registration space (affine)
    mat0a           = seg8t.Affine\M1;            % mat0 for affine output
    mat0r           = R\M1;                       % mat0 for rigid ouput

    % settings for new boundary box of the output images (structure for cat_io_writenii)
    trans.native.Vo = seg8t.image(1); 
    trans.native.Vi = seg8t.image(1);
    trans.affine    = struct('odim',tdim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',seg8t.Affine);  
    trans.rigid     = struct('odim',tdim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);          
  
    job.output.bonemarrow  = struct('native',1,'warped',0,'dartel',3);
    % midline map also for masking masking
    cat_io_writenii(Vo,Ybonemarrow,'',sprintf('bone%d_',job.opts.bmethod), ...
      'bone','uint16',[0,0.001], ... 
      min([1 0 2],[job.output.bonemarrow.native job.output.bonemarrow.warped job.output.bonemarrow.dartel]),trans);
    clear trans
  end
end

 