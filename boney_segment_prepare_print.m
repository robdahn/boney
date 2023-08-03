function [Theader,Tline,Tavg, Cheader, MAfn, matm, mmatm] = ...
  boney_segment_prepare_print(P,job)
%prepare_print. Create table elements for the in-line report and fst print.
%
% [Theader,Tline,Tavg, Cheader, MAfn, matm, mmatm] = ...
%   boney_segment_prepare_print(P,job)
%  
%  P          .. input files
%  job        .. SPM batch structure 
%  Theader    .. header line of the table with all variables names
%  Tline      .. empty subject column 
%  Tavg       .. empty average column
%  Cheader    .. cell header 
%  MAfn       .. cell table with header
%  matm       .. main result table for values
%  mmatm      .. main result table for marks/ratings
% _________________________________________________________________________
%
% Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


% TODO: 
% * final update of table entries and variable names
% * add further parameter in overview


  % Cell structure to define the table header fieldnames and the datatype
  % of the variable and if a seperator (|) should be used. 
% ############# RD202308: We have update the entries and naminges with the finally selected most relevant measures.   
  if 0
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'muscle','fat',   'bone','marrow',   'MBI' , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'  ,   'f'   ,'f'     ,   'f'   , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0    ,   1     ,0       ,   0     , 0    , 0     ]; % 0    ,
  else
    MAfn    = {'Tw','Thbg','Tfat','Tbone','Tres','Tcnr',      'BG',  'CSF','GM',   'Tskull','Thead',  'MED' , 'MBI'  , 'BT' , 'HDT' }; % 'Pll', Tw, Tfat, Thg, Tcnr
    MAffn   = {'s' ,'s'   ,'s'   ,'s'    ,'f'   ,'f'   ,      'f' ,  'f'  ,'f' ,   'f'     ,'f'    ,  'f'   , 'f'    , 'f'  , 'f'   }; % 'e'  ,
    MAfsep  = [0   ,0     ,0     ,0      ,0     ,0     ,      1   ,  0    ,0   ,   1       ,0      ,  1     ,  0     , 0    , 0 ]; % 0    ,
  end
 

  % Creat a string that can be used with the sprintf function to format the
  % different table variables in a nice way.
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',job.opts.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',job.opts.snspace(1)-8);
  Tavg    = sprintf('%%%ds:',job.opts.snspace(1)-1);
  for fi = 1:numel(MAfn)
    Cheader = [Cheader MAfn{fi}]; %#ok<AGROW> 
    AMfni   = strrep( MAfn{fi} ,'_','');
    if MAfsep(fi), sep = ' |'; else, sep = ''; end
    Theader = sprintf(sprintf('%%s%%s%%%ds' ,job.opts.snspace(2)),Theader, sep, AMfni ); 
    switch MAffn{fi}
      case {'s','d'}
        Tline   = sprintf('%s%s%%%d%s', Tline , sep, job.opts.snspace(2), MAffn{fi});
        Tavg    = sprintf('%s%s%%%d%s', Tavg  , sep, job.opts.snspace(2), MAffn{fi});
      otherwise
        Tline   = sprintf('%s%s%%%d.%d%s',Tline , sep, job.opts.snspace(2), job.opts.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
        Tavg    = sprintf('%s%s%%%d.%d%s',Tavg  , sep, job.opts.snspace(2), job.opts.snspace(3) .* (MAffn{fi}=='f'),MAffn{fi});
    end
  end
  % add field for the subject processing time 
  Tline   = sprintf('%s%%4.0fs',Tline);
  

  % main result table
  matm    = num2cell(nan(numel(P),numel(MAfn)));
  mmatm   = 10.5*ones(numel(P),numel(MAfn));
  

  % print title
  if job.opts.verb
    fprintf('\nBone Preprocessing:\n');
    if job.opts.verb 
      methodstr = {'SPMmat8','MedBoneSurf','rMedBoneSurf'};
      reportstr = {'Table','Table + Volumes','Table + Volumes + Surfaces'};
      fprintf('  Method:   %d (%s)\n', job.opts.bmethod, methodstr{job.opts.bmethod+1});
      fprintf('  Report:   %d (%s)\n', job.output.report, reportstr{job.output.report});
      if job.opts.bmethod>0
        fprintf('  Reduce:   %d (%s)\n', job.opts.reduce, reportstr{job.opts.reduce});
      end
      % further parameter ?!
      % ################################
      %  - preprocessing (SEG-CT/SPM/CAT)
      %  - bmethod id + bmethod name 
      %  - report option + name it (only tables, imags, full)
      %  - warnings/information for atypical settings
      % ################################
    end
    fprintf('\n%s\n%s\n',  Theader,repmat('-',size(Theader)));  
  end
end
