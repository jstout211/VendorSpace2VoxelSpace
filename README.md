# AddStandardFids2MOUS

This code moves the fiducials localized on the CTF space MRI to the original T1w.nii.  
This process only works because the data provides the original and ctf space with the same matrix dimensions.

Currently mne-bids-pipeline will error if the ctf space mris are present. <br>

```
cd BIDSDIR; 
mkdir CTFMRIs;
mv sub-*/anat/*CTF_T1w.* CTFMRIs

mne_bids_pipeline/run.py --config=....
```

