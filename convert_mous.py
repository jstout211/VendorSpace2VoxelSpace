#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:14:02 2022

@author: jstout
"""
import os
import os.path as op
import json
import copy
import nibabel as nb
from mne.transforms import apply_trans, invert_transform

topdir = '/fast/OPEN/MOUS_TMP/test_subjs'
os.chdir(topdir)

coordsys_fname='/fast/OPEN/MOUS_TMP/test_subjs/sub-A2002/meg/sub-A2002_coordsystem.json'

	"HeadCoilCoordinateUnits": "cm",
	"HeadCoilCoordinates": {
		"nasion": [10.5092,0,0],
		"left_ear": [-0.014438,7.03839,0],
		"right_ear": [0.014438,-7.03839,0]
	}
    
with open(coordsys_fname) as w:
    coordsys = json.load(w)

def list_nparr2list(nparray):
    return [float(i) for i in nparray]
    
def convert_headcoils2mm(coordsys):
    '''If the fiducials are in cm convert to mm'''    
    coordsys = copy.copy(coordsys)
    #Convert to mm
    if coordsys["HeadCoilCoordinateUnits"].lower()=='cm':
        fids = {key:np.array(value)*10 for key,value in coordsys["HeadCoilCoordinates"].items()}
        fids = {key:list_np2list(value) for key,value in fids.items()}
        coordsys["HeadCoilCoordinates"]=fids
        coordsys["HeadCoilCoordinateUnits"]='mm'
    return coordsys


bids_root='/fast/OPEN/MOUS_TMP/test_subjs'
subject = 'sub-A2002'
intend_for = coordsys['IntendedFor']
intended_for = f'{bids_root}/{subject}/{intend_for}'

anat_t1w = intended_for.replace('space-CTF_','')
affine = nb.load(intended_for).affine

if 'space-CTF' in op.basename(intended_for).split('_'):
    print('True')
    

def make_ras_affine():
    tmp=copy.copy(affine)
    tmp[0,[0,1]]=[-tmp[0,1],-tmp[0,0]]      
    tmp[1,[0,1]]=[-tmp[1,1],-tmp[1,0]]      

        
## The fids are defined in CTF_space
## Need to take ctf affine, invert
## location >> VOX
##
tmp = mne.Transform('ctf_meg', 'mri_voxel', affine)
inv_trans = invert_transform(tmp)
fid_arr=np.stack([np.array(fids['nasion']), 
         np.array(fids['left_ear']),
         np.array(fids['right_ear'])])

tmp_ = np.zeros([3,3])
tmp_[:,0]=-fid_arr[:,1]
tmp_[:,1]=-fid_arr[:,0]
tmp_[:,2]=fid_arr[:,2]
fid_arr=tmp_


fid_vox=apply_trans(inv_trans, fid_arr )
_t1_trans = nb.load(anat_t1w).affine
t1_trans=mne.Transform('mri_voxel', 'mri', _t1_trans)
fid_t1_loc = apply_trans(t1_trans, fid_vox)


# =============================================================================
# Freeview is actually reading out 
# =============================================================================
#CTFmri
ctf_mri_fname=intended_for
ctf_mri=nb.load(ctf_mri_fname)
aff = ctf_mri.affine
mat = ctf_mri.get_fdata()
inv_trans = invert_transform(mne.Transform('ctf_meg','mri_voxel', aff))
fid_vox = apply_trans(inv_trans, fid_arr)

new_fid_vox = np.zeros([3,3])
new_fid_vox[:,0]=fid_vox[:,1]
new_fid_vox[:,1]=fid_vox[:,0]
new_fid_vox[:,2]=fid_vox[:,2]

new_aff = np.zeros([4,4])
new_aff[0,:]=aff[1,:]
new_aff[1,:]=aff[0,:]
new_aff[2,:]=aff[2,:]
new_aff[3,:]=aff[3,:]

#Get the vox
tmp = mne.Transform('ctf_meg', 'mri_voxel', new_aff)
inv_trans = invert_transform(tmp)
fid_arr=np.stack([np.array(fids['nasion']), 
         np.array(fids['left_ear']),
         np.array(fids['right_ear'])])
fid_vox=apply_trans(inv_trans, fid_arr )
_t1_trans = nb.load(anat_t1w).affine
t1_trans=mne.Transform('mri_voxel', 'mri', _t1_trans)
fid_t1_loc = apply_trans(t1_trans, new_fid_vox)


# =============================================================================
# Compare data mats
# =============================================================================
t1 = nb.load(anat_t1w)
t1_ctf = nb.load(ctf_mri_fname)

t1mat = t1.get_fdata()
t1ctfmat = t1_ctf.get_fdata()
t1ctfmat = t1ctfmat[::-1,:,:] 
t1ctfmat = t1ctfmat.transpose(2,0,1)

#T1 index1 = L/R  -- P/A  -- I/S
#T1ctf index1 = AP  -- I/S -- L/R
idx_frame=60
pylab.subplot(2,2,1)
pylab.imshow(t1mat[idx_frame, :, :])
pylab.subplot(2,2,2)
pylab.imshow(t1mat[:, idx_frame, :])
pylab.subplot(2,2,3)
pylab.imshow(t1ctfmat[idx_frame,: ,:])
pylab.subplot(2,2,4)
pylab.imshow(t1ctfmat[:,idx_frame,:])

def convert_ctf2ras(fidval, ctfmat):
    '''Provide the voxel index
    Assumes that both input and output mats are the same size'''
    i0,i1,i2=fidval
    l0,l1,l2 = ctfmat.shape
    i0=l0-i0
    o0 = i2
    o1 = i0
    o2 = i1
    return o0,o1,o2 







    


# =============================================================================
# TESTS
# =============================================================================
def test_unitconv():
    json=	{"HeadCoilCoordinateUnits": "cm",
    	"HeadCoilCoordinates": {
    		"nasion": [10.5092,0,0],
    		"left_ear": [-0.014438,7.03839,0],
    		"right_ear": [0.014438,-7.03839,0]
    	}}
    gtruth =  {'HeadCoilCoordinateUnits': 'mm',
     'HeadCoilCoordinates': {'nasion': [105.092, 0.0, 0.0],
      'left_ear': [-0.14438, 70.3839, 0.0],
      'right_ear': [0.14438, -70.3839, 0.0]}}
    assert convert_headcoils2mm(json)==gtruth
    
def test_ctf2ras():
    nasV = [32,122,87]
    lpaV = [137,93,22]
    rpaV = [130,96,163]
    nasT1V=convert_ctf2ras(nasV, t1ctfmat)
    lpaT1V=convert_ctf2ras(lpaV, t1ctfmat)
    rpaT1V=convert_ctf2ras(rpaV, t1ctfmat)
    
    assert nasT1V==(87,224, 122)
    assert lpaT1V==(22,119,93)
    assert rpaT1V==(163, 126, 96)