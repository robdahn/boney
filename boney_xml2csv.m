function varagout = boney_xml2csv(job)
%boney_xml2csv(job). Wrapper to call cat_io_xml2csv for boney structure. 
% See cat_io_xml2csv.
%
%  varagout = boney_xml2csv(job)
% 
%  
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________
% $Id$

  %job = cat_io_checkin(job,def);

  % define general boney fields that should be avoided even for
  % developers such as help or log entries
  job.avoidfields = [ job.avoidfields; 
    {'help'; 'catlog'; 'error'; 'hardware'; 'P.orgpp';'P.orgff'; }];
  job.dimlim = 10;


  % maybe we need some predefined structures for better named regions etc.  
  job.fieldnames = { ...
    ... filename
    'P.org';
    ... opts
    'opts.pmethod'; 'opts.bmethod';
    'date'; 
    ... SPM image parameter
    'tis.weightingn'; 'tis.highBGn';    'tis.headFatTypen'; 'tis.boneIntTypen';         % basic information about the image
    'tis.bonecortex'; 'tis.bonemarrow'; 'tis.bonedensity';  'tis.bone'; 'tis.head';   
    ... volume measurements
    'tismri.TIV'; 
    'tismri.volr(01)'; 'tismri.volr(02)'; 'tismri.volr(03)'; 'tismri.volr(04)'; 'tismri.volr(05)'; % relative tissue volumes 
    'tismri.volmus'; 'tismri.volfat'; 'tismri.volfatr'; 'tismri.volmusr';
    ... voxel-based measurements (3-occ,4-rpar,5-lpar)
    'vROI.bonecortex(03)';    'vROI.bonemarrow(04)';    'vROI.bonemarrow(05)'; 
    'vROI.bonethickness(03)'; 'vROI.bonethickness(04)'; 'vROI.bonethickness(05)';
    'vROI.headthickness(03)'; 'vROI.headthickness(04)'; 'vROI.headthickness(05)';
    ... surface-based measurements
    'sROI.bonecortex(03)';    'sROI.bonemarrow(04)';    'sROI.bonecortex(05)'; 
    'sROI.bonethickness(03)'; 'sROI.bonethickness(04)'; 'sROI.bonethickness(05)';
    'sROI.headthickness(03)'; 'sROI.headthickness(04)'; 'sROI.headthickness(05)';
    };

 
  % add extra expert/developer fields
  switch job.report
    case {'boney_details','boney_expertx','boney_developer'}
      job.fieldnames = [job.fieldnames; {
        'tismri.vol(01)';  'tismri.vol(02)';  'tismri.vol(03)';  'tismri.vol(04)';  'tismri.vol(05)';  % absolute tissue volumes
        'vROI.bonemarrow';'vROI.bonecortex';'vROI.bonethickness';'vROI.headthickness';'vROI.head';
        'sROI.bonecortex';'sROI.bonecortex';'sROI.bonethickness';'sROI.headthickness';'sROI.head';
        } ]; 
  end
    
  switch job.report
    case {'boney_expertx','boney_developer'}
      job.fieldnames = [job.fieldnames; {
        'opts.bias';'opts.ctpm'; 'opts.refine'; 'opts.reduce';
        } ]; 
  end
  switch job.report
    case {'boney_developer'}
      job.fieldnames = [job.fieldnames; {
        'tis.clsQC';
        } ]; 
  end

  varagout{1} = cat_io_xml2csv(job);

end
