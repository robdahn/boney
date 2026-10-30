function varargout = boney_segment_QC(action,varargin)
%boney_segment_QC. Quick quality control flags of the boney processing.
% Simple tests of the SPM preprocessing (seg8 values, affine registration)
% and of the boney measures to detect typical problems such as failed
% segmentations, skull-stripped or defaced images, or implausible values.
% Serious problems are reported as warnings with a suggestion.
%
%  tis = boney_segment_QC('seg8',tis,seg8t)   % SPM-based tests
%  qc  = boney_segment_QC('final',out)        % all tests (incl. MRI-based)
%  fail = boney_segment_QC('brainseg',volr)   % test of GM/WM/CSF of TIV
%
%  qc         .. QC structure
%   .flags    .. short names of the detected problems (cellstr)
%   .fail     .. likely failed processing (severe problem)
%   .segfail  .. likely failed segmentation that may be solved by another
%                TPM (used for the automatic TPM fallback in boney_segment)
%   .val      .. values used for the tests
%   .tpmsel   .. selected TPM in case of automatic selection
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________

  switch action
    case 'seg8'
      varargout{1} = qc_seg8(varargin{:});
    case 'final'
      varargout{1} = qc_final(varargin{:});
    case 'brainseg'
      varargout{1} = qc_brainseg(varargin{:});
    otherwise
      error('boney:QC:action','Unknown action "%s".',action);
  end
end
%=======================================================================
function tis = qc_seg8(tis,seg8t)
%qc_seg8. Tests based on the SPM seg8 values.

  qc.flags = {}; qc.fail = false;

  % affine registration: isotropic scaling (head size relative to the TPM)
  % and rotation - large values suggest a failed affine registration
  if isfield(seg8t,'Affine') && all(isfinite(seg8t.Affine(:)))
    P = spm_imatrix(seg8t.Affine);
    qc.val.affzoom = abs(det(seg8t.Affine(1:3,1:3)))^(1/3);
    qc.val.affrot  = max(abs(P(4:6))) * 180/pi;
    qc.val.affaniso = max(abs(P(7:9))) / min(abs(P(7:9)));
    if qc.val.affzoom < 0.75 || qc.val.affzoom > 1.4 || qc.val.affrot > 45 || qc.val.affaniso > 1.4
      qc.flags{end+1} = 'affine'; qc.fail = true;
    end
  end

  % image contrast: the bone/head measures are optimized for T1w images
  if tis.weighting ~= 1
    qc.flags{end+1} = sprintf('contrast%s',tis.weightingn);
  end

  % high noise (WM SD in percent of the WM intensity)
  qc.val.CNR = tis.seg8CNR;
  if tis.seg8CNR > 15
    qc.flags{end+1} = 'noise';
  end

  tis.qc = qc;
end
%=======================================================================
function fail = qc_brainseg(volr)
%qc_brainseg. Relative brain tissue volumes (GM, WM, CSF of the TIV) - a
% collapsed or overestimated tissue class indicates a failed segmentation
% (wide limits as elderly subjects have much CSF and the BlaiottaTPM
% results in less WM than the SPM TPM).
  fail = volr(1) < 0.20 || volr(1) > 0.70 || volr(2) < 0.15 || volr(2) > 0.55 || volr(3) > 0.50;
end
%=======================================================================
function qc = qc_final(out)
%qc_final. Combine the SPM-based tests with tests of the MRI-based values.

  if isfield(out,'tis') && isfield(out.tis,'qc')
    qc = out.tis.qc;
  else
    qc.flags = {}; qc.fail = false;
  end

  if isfield(out,'tismri') && isfield(out.tismri,'TIV')
    tm = out.tismri;

    % TIV - typical range of adults and children
    qc.val.TIV = tm.TIV;
    if tm.TIV < 700 || tm.TIV > 2300
      qc.flags{end+1} = 'TIV'; qc.fail = true;
    end

    % relative brain tissue volumes (GM, WM, CSF)
    qc.val.volr = tm.volr(1:3);
    if qc_brainseg(tm.volr(1:3))
      qc.flags{end+1} = 'brainseg'; qc.fail = true;
    end

    if numel(tm.volr) > 4
      % skull-stripped images: almost no bone and head
      qc.val.volrBoneHead = tm.volr(4) + tm.volr(5);
      if tm.volr(4) + tm.volr(5) < 0.1
        qc.flags{end+1} = 'skullstripped'; qc.fail = true;
      end

      % implausible bone class: the bone dominates the extracerebral 
      % tissues and replaces the CSF, i.e., CSF and head tissues were 
      % labeled as bone (e.g. children or non-T1w low-resolution scans 
      % with the SPM TPM). The bone is typically <30% of the bone and head
      % volume (max. 0.45 for CSF<0.2 in the development, UKB, and NKI 
      % data), whereas the failed cases showed >0.45 with CSF<0.2.
      qc.val.boneHeadRatio = tm.volr(4) / max(eps,tm.volr(4) + tm.volr(5));
      if qc.val.boneHeadRatio > 0.45 && tm.volr(3) < 0.2
        qc.flags{end+1} = 'boneclass'; qc.fail = true;
      end
    end
  end

  % defacing or masked background (ratio of the outer head boundary with 
  % zero intensity)
  if isfield(out,'tismri') && isfield(out.tismri,'zeroBGratio')
    qc.val.zeroBGratio = out.tismri.zeroBGratio; 
    if out.tismri.zeroBGratio > 0.8
      qc.flags{end+1} = 'maskedBG'; 
    elseif out.tismri.zeroBGratio > 0.03
      qc.flags{end+1} = 'defaced'; 
    end
  end

  % bone measures (the global bone thickness of successful cases was 
  % below 11.5 mm in the development, UKB, and NKI data)
  if isfield(out,'vROI') && isfield(out.vROI,'bonethickness')
    qc.val.bonethickness = out.vROI(1).bonethickness(1);
    if isnan(out.vROI(1).bonethickness(1)) || out.vROI(1).bonethickness(1) < 2 || out.vROI(1).bonethickness(1) > 12
      qc.flags{end+1} = 'bonethickness'; qc.fail = true;
    end
  end
  if isfield(out,'main') && isfield(out.main,'sBMDH') && isnan(out.main.sBMDH)
    qc.flags{end+1} = 'noBMD';
  end

  % automatic TPM selection
  if isfield(out,'tpmsel') && isfield(out.tpmsel,'ctpm')
    qc.tpmsel = out.tpmsel.ctpm;
  else
    qc.tpmsel = nan;
  end

  qc.flags   = unique(qc.flags,'stable');
  qc.segfail = any(ismember(qc.flags,{'TIV','brainseg','boneclass','bonethickness'}));
  qc.help  = ['Quick QC flags: affine (implausible affine registration), contrast* (not T1w), ' ...
    'noise (WM SD>15%), TIV (TIV<700 or >2300 cm3), brainseg (implausible GM/WM/CSF volumes), ' ...
    'skullstripped (no bone/head), boneclass (bone >45% of bone and head with CSF <20% of TIV), ' ...
    'defaced/maskedBG (zero intensities at the head boundary), bonethickness (<2 or >12 mm), ' ...
    'noBMD (no surface BMD). fail marks likely failed processing, segfail a failed segmentation ' ...
    'that is processed again with another TPM in case of the automatic TPM selection.'];

  % warning with suggestion
  if qc.fail
    cat_io_addwarning('boney:QC:fail', sprintf(['Likely failed processing (QC flags: %s). ' ...
      'Check the report and try another TPM (or the automatic TPM selection).'], strjoin(qc.flags,', ')), 2, [1 1], {}, 0, 0);
  end
end
