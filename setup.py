import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VendorSpace2VoxelSpace", 
    version="0.1",
    author="Jeff Stout",
    author_email="stoutjd@nih.gov",
    description="Convert non-RAS defined fiducials into voxel based fids",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jstout211/VendorSpace2VoxelSpace",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: UNLICENSE",
        "Operating System :: Linux/Unix",
    ],
    install_requires=['mne', 'numpy', 'scipy', 'pandas', 'nibabel', 'pytest', 'joblib', 'mne_bids'],
    #scripts=['process_meg.py', 
    #         'process_anatomical.py'],
)
