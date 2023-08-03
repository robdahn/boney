# Boney
This is the Boney toolbox, an extension to [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and its [CAT12](http://www.neuro.uni-jena.de/cat) toolbox to support additional measurements to quantify the bone and head structure.  It is developed by *Polona Kalc* and *Robert Dahnke* and free but copyright software, distributed under the terms of the <em>[GNU General Public License](http://www.gnu.org/licenses/gpl-2.0.html)</em> as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

If you find any bug, please report them to <polona.kalc@med.uni-jena.de> or <robert.dahnke@uni-jena.de>.


## Introduction 
Why some people believe that bones are more interesting than a naked women ...
... *bone mineral density* (**BMD**) 
... *body mass index* (**BMI**) 


## Quick method
To extract bone parameters from MRI data the SPM tissue segmetnation is estimated and refined to extract the different bone and head intensity and thickness measures that can be used to approximate the BMD or the BMI of a set of subjects. 


## Quick start
Install Matlab, [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and its [CAT12](http://www.neuro.uni-jena.de/cat).
Download the zip-file and unpack it into your SPM12 toolbox directory. 

Run SPM and open the toolbox:

<code>spm fmri; spm_boney;</code>

Open the bone processing batch: 
![Image of the Boney menu and the bone-processing-batch](/images/Boney_toolbox.png "Shown is the bone-processing batch that can be used to extract processed bone measures into a CSV table.")

and select the structural images that should be processed, adopt some parameters (see paragraph Parameters and SPM batch help) and start the processing. 

After the processing has finished, you can use the XML2CSV batch to extract the bone values of relevant subjects within one table for further analyses.
![Image of the Boney menu and the XML2CSV-batch](/images/boney_software_XML2CSV.jpg "Shown is the XML2CSV batch that can be used to extract processed bone measures into a CSV table.")


## Parameters
You can select different predefined SPM and CAT segmentation routines with different predefined settings. 
SPM is faster but may fail in some cases, often by inital affine registration problems.

In addition you can select between different complex bone processes routines - SPMmat, volume, surface.  
...

Furthermore, you can specify whenever output files are writen, such as the short bone-report (as a JPG), or the processed NIFTI volumes or GIFTI surfaces that are also shown in the bone-report.


## Results Validation 
We used the UKB to ... 


## References



