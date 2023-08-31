function Pout = boney_segment_cleanup(Pout,out,job,i)
%boney_segment_cleanup. Update Pout structure and clean up files  
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

  if job.output.writeseg == 1
    for ci = 1:5
      Pout.cls{ci}{i}    = out(i).P.cls{ci}; 
    end
  elseif 0
    % delete preprocessing segments
    for ci = 1:5
      if exist( out(i).P.cls{ci} , 'file')
        delete( out(i).P.cls{ci} );
      end
    end
  end


  % boney volume files
  if job.output.writevol
    Pout.rbone_affine{i} = out(i).P.bone_affine; 
    Pout.rbone{i}        = out(i).P.bone;
  else
    % delete files
  end


  % boney surface files
  if job.output.writesurf
    Pout.bcentral{i}     = out(i).P.central;  
    Pout.bcortex{i}      = out(i).P.cortex;
    Pout.bmarrow{i}      = out(i).P.marrow;
    Pout.bthickness{i}   = out(i).P.thick; 
    Pout.hthickness{i}   = out(i).P.headthick; 
  elseif job.opts.bmethod == 2 % surface processing
    FD = {'central','cortex','marrow','thick','headthick'}; 
    for fdi = 1:numel(FD)
      if exist( out(i).P.(FD{fdi}) , 'file')
        delete( out(i).P.(FD{fdi}) );
      end
    end
  end

end
