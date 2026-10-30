function Pout = boney_segment_cleanup(Pout,out,job,i,ppcreated)
%boney_segment_cleanup. Update Pout structure and clean up files  
%
%  Pout = boney_segment_cleanup(Pout,out,job,i,ppcreated)
%
%  ppcreated .. the SPM/CAT preprocessing files were created by boney in 
%               this run and can be removed (files given by the user or 
%               from previous runs are never removed)
% _________________________________________________________________________
%
% Robert Dahnke & Polona Kalc
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% _________________________________________________________________________


  % xml file and report
  Pout.xml{i}            = out(i).P.xml;
  Pout.reports{i}        = out(i).P.report; 

  if ~exist('ppcreated','var'), ppcreated = false; end
  if job.output.writeseg == 1 || ~ppcreated
    for ci = 1:numel(out(i).P.cls) % 6 classes for BlaiottaTPM
      Pout.cls{ci}{i}    = out(i).P.cls{ci};
    end
  else 
    % delete preprocessing segments (in case of writeseg == 2 they are 
    % stored in the boneyPPmat file)
    for ci = 1:numel(out(i).P.cls) % 6 classes for BlaiottaTPM
      if exist( out(i).P.cls{ci} , 'file')
        delete( out(i).P.cls{ci} );
      end
    end
    % delete bias corrected
    if exist( out(i).P.bc , 'file')
      delete( out(i).P.bc );
    end
    % delete SPM mat
    if exist( out(i).P.seg8 , 'file') && strcmp(out(i).P.seg8(end-3:end),'.mat')
      delete( out(i).P.seg8 );
    end
  end


  % boney volume files
  % (writevol: 1-native, 2-affine, 3-warped, 4-affine+warped, 5-all)
  if any(job.output.writevol == [1 5])
    Pout.bone{i}         = out(i).P.bone;
  end
  if any(job.output.writevol == [2 4 5])
    Pout.rbone_affine{i} = out(i).P.rbone_affine; 
  end


  % boney surface files
  if job.output.writesurf && job.opts.bmethod == 2
    Pout.bcentral{i}     = out(i).P.central;  
    Pout.bcortex{i}      = out(i).P.cortex;
    Pout.bmarrow{i}      = out(i).P.marrow;
    Pout.bthickness{i}   = out(i).P.thick; 
    Pout.hthickness{i}   = out(i).P.headthick; 
  elseif job.opts.bmethod == 2 % surface processing
    FD = {'central','cortex','marrow','thick','headthick'}; 
    for fdi = 1:numel(FD)
      if isfield(out(i).P,FD{fdi}) && exist( out(i).P.(FD{fdi}) , 'file')
        delete( out(i).P.(FD{fdi}) );
      end
    end
  end

end
