import os
from Dicom4D import Dicom4D

def propagate(fndcm, outdir = "", tag = "", seg_ref = "", framenums = "", fref = ""):

    """
    INPUTS:
    - fndcm: Cartesian DICOM filename
    - outdir: Name of existing directory for output files
    - tag: Study Identifier (e.g., 'bav08_root')
    - seg_ref: Filename of reference segmentation (nii)
    - framenums: List of frame numbers to segment
    - fref: Reference segmentation frame

    Use full paths for all filenames.

    """

    """
    Process Inputs:
    - Load cartesian dicom image
    - Process reference segmentation
        - Dialate the reference segmentation
        - Create vtk mesh based on dialted reference segmentation
    - Extract and dialate 3D frame from the 4D image

    """

    # Use the  reader to create a Dicom4D Image
    CartesianDicom = Dicom4D(fndcm)
    # Export 3D Frames
    # CartesianDicom.Export4D(os.path.join(outdir, 'img4D.nii.gz'))
    for i in framenums:
        fnImg = f'{outdir}/img{i:02d}_{tag}.nii.gz'
        fnImgRs = f'{outdir}/img{i:02d}_{tag}_srs.nii.gz'
        CartesianDicom.ExportFrame(i, fnImg)
        
        syscmd = '/usr/local/bin/c3d ' + fnImg + ' -smooth 1mm -resample 50% \-o ' + fnImgRs
        print(syscmd)
        os.system(syscmd)


    # Propagate Forward
    


    # Propagate Backward



    # Propagate in Full Resolution



