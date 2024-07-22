Pt1 = '/Volumes/OneTouch/MRData/202301_MRIbones/ukb_bones/VBM_women_1000_headBMD/spmT_0002.nii';
Pt2 = '/Volumes/OneTouch/MRData/202301_MRIbones/ukb_bones/VBM_women_1000_BMDH_WM/spmT_0002.nii';

Vt1 = spm_vol(Pt1);
Vt2 = spm_vol(Pt2);

Yt1 = spm_read_vols(Vt1); 
Yt2 = spm_read_vols(Vt2); 

%%
corr(Yt1(Yt1(:)~=0),Yt2(Yt2(:)~=0))
