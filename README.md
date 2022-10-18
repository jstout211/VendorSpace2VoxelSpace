# VendorSpace2VoxelSpace
Convert fiducials in fieldtrip's MEG vendor space bids output to MRI voxel space.


# AddStandardFids2MOUS
This code moves the fiducials localized on the CTF space MRI json in cm to the original T1w.nii json in voxel coordinates.  
This process only works because the data provides the original and ctf space with the same matrix dimensions.

```
convert_mous.py MOUSBIDSDIR
```
Verify that all processes complete successfully before doing the next step. <br>
Currently mne-bids-pipeline will error if the ctf space mris are present. <br>

```
cd MOUSBIDSDIR; 
mkdir CTFMRIs;
mv sub-*/anat/*CTF_T1w.* CTFMRIs

mne_bids_pipeline/run.py --config=....
```

