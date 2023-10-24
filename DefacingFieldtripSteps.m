clear 

addpath(genpath('~\fieldtrip-20220202'))

Patients={'sub-Pt01'};
List={'test'};
%%
clear ct mri mri_anon ct_anon

%%
Direc='M:\XXXXXX\DepthPaper\sub-Pt01_MicrosEEG\ses-preimp\anat\';
FiName=[Direc,'\sub-Pt01_MicrosEEG_ses-postimp_acq-MRI_run-01_T1w.nii'];
mri = ft_read_mri(FiName);
cfg = [];
mri_anon = ft_defacevolume(cfg, mri);
%
ft_write_mri(FiName, mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');

%%
Direc='M:\XXXXXX\DepthPaper\sub-Pt01_MicrosEEG\ses-postimp\anat\';
FiName=[Direc,'\sub-Pt01_MicrosEEG_ses-postimp_acq-CT_run-01_T1w.nii'];
mri = ft_read_mri(FiName);
cfg = [];
mri_anon = ft_defacevolume(cfg, mri);
%
ft_write_mri(FiName, mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');

