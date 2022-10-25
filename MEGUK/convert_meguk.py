#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:14:02 2022

@author: jstout
"""

import logging
logging.basicConfig(
                level=logging.INFO,
                force=True,
                format='%(asctime)s %(message)s',
                handlers=[logging.FileHandler("meguk_fid_conv.log")])

import os
import os.path as op

import json
import copy
import nibabel as nb
import numpy as np
import glob
from mne.transforms import apply_trans, invert_transform
import mne
import mne_bids
from mne_bids import BIDSPath

# def get_fiducials():
    # subjid=subjid.copy()
    # if subjid[0:4]=='sub-':
    #     subjid=subjid[4:]
        
    # if subjid[0:3]=='gla':
    #     #Dataset is in 
    # if subjid[0:3]=='cdf':
    #     coordsys='test'
    # return    



def list_nparr2list(nparray):
    return [float(i) for i in nparray]
    
def convert_headcoils2mm(coordsys):
    '''If the fiducials are in cm convert to mm'''    
    coordsys = copy.copy(coordsys)
    #Convert to mm
    if coordsys["HeadCoilCoordinateUnits"].lower()=='cm':
        fids = {key:np.array(value)*10 for key,value in coordsys["HeadCoilCoordinates"].items()}
        fids = {key:list_nparr2list(value) for key,value in fids.items()}
        coordsys["HeadCoilCoordinates"]=fids
        coordsys["HeadCoilCoordinateUnits"]='mm'
    return coordsys

# def convert_ctf2t1(fidval, ctfmat):
#     '''Provide the voxel index
#     Assumes that both input and output mats are the same size'''
#     i0,i1,i2=fidval
#     l0,l1,l2 = ctfmat.squeeze().shape  #Squeeze fixes singleton dimensions
#     i0=l0-i0
#     o0 = i1
#     o1 = i0
#     o2 = i2
#     return o0,o1,o2 

def write_anat_json(anat_json=None, 
                    fids=None,
                    overwrite=False):
    '''Write the AnatomicalLandmark keys from the meg json to the mri json.
    Takes the dictionary pulled from the MEG json - (use get_anat_json)'''
    with open(anat_json) as w:
        json_in = json.load(w)
        
    #Check to see if the AnatomicalLandmarks are already present 
    anat_keys = {x:y for (x,y) in json_in.items() if 'AnatomicalLandmark' in x}
    if overwrite==False:
        assert anat_keys == {}

    json_in['AnatomicalLandmarkCoordinates']=fids  
    
    #Write json
    with open(anat_json, 'w') as w:
        json.dump(json_in, w, indent='    ') 


def convert_single_subject(bids_root=None,
                           subject=None,
                           overwrite=True
                           ):
    
    coordsys_fname = f'{bids_root}/{subject}/meg/{subject}_coordsystem.json'
    with open(coordsys_fname) as w:
        coordsys = json.load(w)
        
    intend_for = coordsys['IntendedFor']
    intended_for = f'{bids_root}/{subject}/{intend_for}'
    
    if 'space-CTF' not in op.basename(intended_for).split('_'):
        raise('space-CTF is not in the intended for label')
    
    anat_t1w = intended_for.replace('space-CTF_','')
    if os.path.exists(anat_t1w.replace('.nii','.json')):
        anat_t1w_json = anat_t1w.replace('.nii','.json')
    elif os.path.exists(anat_t1w.replace('.nii.gz','.json')):
        anat_t1w_json = anat_t1w.replace('.nii.gz','.json')
    fids = convert_headcoils2mm(coordsys)['HeadCoilCoordinates']
    
    #Get the FIDS
    fid_arr=np.stack([np.array(fids['nasion']), 
             np.array(fids['left_ear']),
             np.array(fids['right_ear'])])
    
    #CTFmri
    ctf_mri_fname=intended_for
    ctf_mri=nb.load(ctf_mri_fname)
    aff = ctf_mri.affine
    mat = ctf_mri.get_fdata()
    inv_trans = invert_transform(mne.Transform('ctf_meg','mri_voxel', aff))
    fid_vox = apply_trans(inv_trans, fid_arr)
    
    nasT1V=convert_ctf2t1(fid_vox[0], mat)
    lpaT1V=convert_ctf2t1(fid_vox[1], mat)
    rpaT1V=convert_ctf2t1(fid_vox[2], mat)
    
    t1_vox_fids=dict(NAS=nasT1V,
                     LPA=lpaT1V,
                     RPA=rpaT1V)
    write_anat_json(anat_json=anat_t1w_json, 
                    fids=t1_vox_fids,
                    overwrite=overwrite)
    logging.info(f'Success: {subject}')
    
    
def convert_mous_project(bids_dir=None):
    '''Loop over all subjects and add the Voxel index Fiducials to the 
    T1w.json file in the anatomy folder
    '''
    bids_subjs = glob.glob(op.join(bids_dir, 'sub-*'))
    bids_subjs = [op.basename(i) for i in bids_subjs] 
    
    for subjid in bids_subjs:
        print(f'Running: {subjid}')
        logging.info(f'Running: {subjid}')
        try:
            convert_single_subject(bids_dir, subject=subjid, overwrite=True)
        except BaseException as e:
            print(f'Error with {subjid}:\n {e}')
            logging.Exception(f'Error with {subjid}:\n {e}')
 

# =============================================================================
#         
# =============================================================================
#%%

def convert_ctf2t1(fidval, ctfmat):
    '''Provide the voxel index
    Assumes that both input and output mats are the same size'''
    i0,i1,i2=fidval
    l0,l1,l2 = ctfmat.squeeze().shape  #Squeeze fixes singleton dimensions
    i0=l0-i0
    o0 = i0
    o1 = i1
    o2 = i2
    return o0,o1,o2 


class fids_space_to_vox():
    '''Convert vendor space to voxel space'''
    def __init__(self, subject=None, bids_root=None, session=1):
        self.subject=subject
        self.bids_root=bids_root
        self.bids_path=BIDSPath(root=bids_root, 
                                subject=subject, 
                                session=session)
        self.fids=None
        self.write_anat_json=write_anat_json
        self.get_coordsys_fname()
        self.get_fids()
        self.get_fids_voxel_ctf()
        self.anat_json_fname = self.bids_path.copy().update(datatype='anat',
                                                            suffix='T1w',
                                                            extension='.json')
        self.anat_json_fname = str(self.anat_json_fname.fpath)
        self.write_anat_json(anat_json=self.anat_json_fname, 
                    fids=self.fids_T1w,
                    overwrite=True)

    def get_coordsys_fname(self):
        _bidspath=self.bids_path.copy()
        if _bidspath.update(datatype='meg', suffix='coordsystem', extension='tsv').fpath:
            self.coordsys_json=_bidspath.update(datatype='meg', suffix='coordsystem', extension='json').fpath
    
    def get_fids(self):    
        with open(self.coordsys_json) as w:
            coordsys = json.load(w)
        self.coordsys = coordsys
            
        intend_for = coordsys['IntendedFor']
        intended_for = f'{self.bids_root}/sub-{self.subject}/{intend_for}'
        self.intended_for = intended_for
        
        if 'space-CTF' not in op.basename(intended_for).split('_'):
            raise('space-CTF is not in the intended for label')
        
        anat_t1w = intended_for.replace('space-CTF_','')
        self.anat_t1w = anat_t1w
        if os.path.exists(anat_t1w.replace('.nii','.json')):
            anat_t1w_json = anat_t1w.replace('.nii','.json')
        elif os.path.exists(anat_t1w.replace('.nii.gz','.json')):
            anat_t1w_json = anat_t1w.replace('.nii.gz','.json')
        self.fids = convert_headcoils2mm(coordsys)['HeadCoilCoordinates']
        
    def get_fids_voxel_ctf(self):
        fids=self.fids
        fid_arr=np.stack([np.array(fids['nasion']), 
                  np.array(fids['left_ear']),
                  np.array(fids['right_ear'])])
        
        #CTFmri
        ctf_mri_fname=self.intended_for
        ctf_mri=nb.load(ctf_mri_fname)
        aff = ctf_mri.affine
        self.ctf_aff = ctf_mri.affine
        mat = ctf_mri.get_fdata()
        self.ctf_mat = mat.squeeze()
        inv_trans = invert_transform(mne.Transform('ctf_meg','mri_voxel', aff))
        self.ctf_mrivox_trans = inv_trans
        fid_vox = apply_trans(inv_trans, fid_arr)
        
        nas_val = convert_ctf2t1(fid_vox[0], mat)
        lpa_val = convert_ctf2t1(fid_vox[1], mat)
        rpa_val = convert_ctf2t1(fid_vox[2], mat)
        
        
        self.fids_vox_ctf = fid_vox
        self.fids_T1w=dict(NAS=[i for i in nas_val], 
                          LPA=[i for i in lpa_val],
                          RPA=[i for i in rpa_val])
        
    def plot_alignment(self):
        raw_bids_path = self.bids_path.copy().update(datatype='meg',
                                                 task='resteyesclosed',
                                                 run=None,
                                                 session=None, 
                                                 extension='.ds')
        raw = mne_bids.read_raw_bids(raw_bids_path).pick_types(meg=True)
        anat_dir=os.path.dirname(self.anat_json_fname)
        tmp_dir=os.path.join(anat_dir, 'TMP')
        import shutil
        try:
            if not os.path.exists(tmp_dir): os.mkdir(tmp_dir)
            fnames_mv = glob.glob(os.path.join(anat_dir, '*space-*'))
            for i in fnames_mv:
                shutil.move(i, tmp_dir)
        except:
            print('failed')
        self.t1w_bids_path = raw_bids_path.copy().update(datatype='anat', task=None, suffix='T1w', extension='.nii.gz', 
                                                             space=None)
        trans = mne_bids.get_head_mri_trans(raw_bids_path, 
                                                t1_bids_path=self.t1w_bids_path,
                                                fs_subject='sub-'+self.bids_path.subject) 
        mne.viz.plot_alignment(raw.info,
                                    trans=trans,
                                    subject='sub-'+self.bids_path.subject, dig=True)
        fnames_mv = glob.glob(os.path.join(tmp_dir, '*'))
        if fnames_mv != []:
            for i in fnames_mv:
                shutil.move(i, anat_dir)
        
    def plot_mriview(self):
        import pylab
        import nibabel as nb
        self.t1mrimat = nb.load(self.anat_t1w).get_fdata()
        
        pylab.subplot(2,3,1)
        pylab.imshow(self.ctf_mat[int(self.fids_vox_ctf[0,0]), :, :])
        pylab.subplot(2,3,2)
        pylab.imshow(self.ctf_mat[:, int(self.fids_vox_ctf[0,1]), :])
        pylab.subplot(2,3,3)
        pylab.imshow(self.ctf_mat[:, :, int(self.fids_vox_ctf[0,2])])
        
        pylab.subplot(2,3,4)
        pylab.imshow(self.t1mrimat[int(self.fids_T1w['NAS'][0]), :, :])
        pylab.subplot(2,3,5)
        pylab.imshow(self.t1mrimat[:, int(self.fids_T1w['NAS'][1]), :])
        pylab.subplot(2,3,6)
        pylab.imshow(self.t1mrimat[:, :, int(self.fids_T1w['NAS'][2])])
                               
def test_meg_uk(subjid):
    os.chdir('/fast/EnigmaTesting')
    return fids_space_to_vox(subject=subjid, bids_root='/fast/EnigmaTesting/bids',session=None) 

# for i in ['sub-cdf040', 'sub-cdf043', 'sub-cdf046', 'sub-cdf044', 'sub-cdf045', 'sub-cdf049', 'sub-cdf041', 'sub-cdf047', 'sub-cdf048', 'sub-cdf042']:
#     subjid = i.split('-')[1]
#     test_meg_uk(subjid)


# tmp = test_meg_uk()
# tmp.plot_mriview()
# tmp.plot_alignment()


#%%
os.environ['SUBJECTS_DIR']='/fast/EnigmaTesting/bids/derivatives/freesurfer/subjects'

tmp_40=test_meg_uk('cdf040')
tmp_43=test_meg_uk('cdf043')
tmp_46=test_meg_uk('cdf046')
tmp_44=test_meg_uk('cdf044')
tmp_45=test_meg_uk('cdf045')
tmp_49=test_meg_uk('cdf049')
tmp_41=test_meg_uk('cdf041')
tmp_47=test_meg_uk('cdf047')
tmp_48=test_meg_uk('cdf048')
tmp_42=test_meg_uk('cdf042')


# =============================================================================
# 
# =============================================================================


if __name__=='__main__':
    import sys
    bids_dir = sys.argv[1]
    convert_mous_project(bids_dir=bids_dir)

# =============================================================================
# TESTS MEGUK
# =============================================================================





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
    
def test_ctf2t1():
    nasV = [32,122,87]
    lpaV = [137,93,22]
    rpaV = [130,96,163]
    nasT1V=convert_ctf2t1(nasV, t1ctfmat)
    lpaT1V=convert_ctf2t1(lpaV, t1ctfmat)
    rpaT1V=convert_ctf2t1(rpaV, t1ctfmat)
    
    assert nasT1V==(87,224, 122)
    assert lpaT1V==(22,119,93)
    assert rpaT1V==(163, 126, 96)