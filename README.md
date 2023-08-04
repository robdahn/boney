# Boney
This is the Boney toolbox, an extension to [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and its [CAT12](http://www.neuro.uni-jena.de/cat) toolbox, supporting the extraction of additional measurements related to the bone and head structure. It is developed by *Polona Kalc* and *Robert Dahnke* and is a free but copyright software, distributed under the terms of the <em>[GNU General Public License](http://www.gnu.org/licenses/gpl-2.0.html)</em> as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

If you find any bugs, please report them to <polona.kalc@med.uni-jena.de> or <robert.dahnke@uni-jena.de>.


## Introduction 
As Ingmar Bergman posited in The Seventh Seal: *'A skull is more interesting than a naked woman.'* Yet the neuroimaging community has continued to strip the skull and discard the relevant information hidden in the layers of the head tissue. 
The cranium is the most proximate bone tissue to the brain, and recent studies suggest an intensive bidirectional communication of bones and the brain, and the implication of the bone metabolism in mood, cognition, and energy metabolism ([1]Guntur and Rosen 2012; Khrimian et al. 2017; Nakamura, Imaoka, and Takeda 2021; Obri et al. 2018; Oury et al. 2013; Rousseaud et al. 2016). Furthermore, low bone mineral density (BMD) has been associated to an increased risk of Alzheimer's disease (refs).
However, the bone mineral density measure is typically not available in open-access brain-imaging databases, such as IXI, ADNI, AIBL, OASIS, etc. We therefore present proxy measures for head BMD estimation and head skin-fold thickness that can be used as an additional information ........


## Quick method
To extract bone parameters from MRI data, we use a refined SPM tissue segmentation procedure and derive different cranial intensity and thickness measures that can be used to approximate the head BMD measure.
An estimation of the head's skin-fold (subcutaneous fat thickness) is also available and can be used as an additional information to the typically available BMI measure.


## Quick start
Install Matlab/Octave, [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and its [CAT12](http://www.neuro.uni-jena.de/cat) toolbox.
Download the zip-file and unpack it into your SPM12 toolbox directory. 

Run SPM and open the toolbox:

<code>spm fmri; spm_boney;</code>

Open the bone processing batch: 
![Image of the Boney menu and the bone-processing-batch](/images/Boney_toolbox.png "Shown is the bone-processing batch that can be used to extract processed bone measures into a CSV table.")

and select the structural images that should be processed. Adapt relevant parameters (see paragraph Parameters and SPM batch help) and start the processing. Voilà! 

After the processing has finished, you can use the XML2CSV batch to extract the estimated (bone) values for further analyses.



## Parameters
You can select different predefined SPM and CAT segmentation routines with various predefined settings. 
SPM's routine is faster but may fail in some cases, often because of the inital affine registration problems. 

In addition, you can select between different advanced processing routines for the bone measures (e.g., SPMmat, volume, surface).  
We would advise the use of the refined measure, which is robust to the exclusion of the parts within the diploë.

Furthermore, you can specify which output files to write, e.g., the short bone-report (as a JPG), or the processed NIFTI volumes or GIFTI surfaces that are also shown in the bone-report.

## References
You can find out more about the bone and skin-fold thickness measures in the paper ...


