import os
import vtk
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
        - Dilate the reference segmentation
        - Create vtk mesh based on dialted reference segmentation
    - Extract and dialate 3D frame from the 4D image

    """

    # directory storing temporary files
    tmpdir = os.path.join(outdir, "tmp")

    # Use the  reader to create a Dicom4D Image
    CartesianDicom = Dicom4D(fndcm)

    # Process reference segmentation
    # - Dilate the reference segmentation (mask)
    fn_mask_ref_srs = os.path.join(tmpdir, f'mask_{fref:02d}_{tag}_srs.nii.gz')
    cmd = f"/usr/local/bin/c3d -int 0 {seg_ref} -threshold 1 inf 1 0 \
        -dilate 1 10x10x10vox -resample 50% -o {fn_mask_ref_srs}"
    print("Dilating reference segmentation...")
    print(cmd)
    os.system(cmd)

    # - Create vtk mesh
    fn_mask_ref_vtk = os.path.join(tmpdir, f'mask_{fref:02d}_{tag}.vtk')
    cmd = f"vtklevelset {seg_ref} {fn_mask_ref_vtk} 1"
    print("Creating mask mesh...")
    print(cmd)
    os.system(cmd)

    # - Create vtk mesh from the dilated mask
    fn_mask_ref_srs_vtk = os.path.join(tmpdir, f'mask_{fref:02d}_{tag}_srs.vtk')
    cmd = f"vtklevelset {fn_mask_ref_srs} {fn_mask_ref_srs_vtk} 1"
    print("Creating dilated mask mesh...")
    print(cmd)
    os.system(cmd)



    # Export 3D Frames
    # CartesianDicom.Export4D(os.path.join(outdir, 'img4D.nii.gz'))
    for i in framenums:
        fnImg = f'{outdir}/img{i:02d}_{tag}.nii.gz'
        fnImgRs = f'{outdir}/img{i:02d}_{tag}_srs.nii.gz'
        CartesianDicom.ExportFrame(i, fnImg)
        
        cmd = '/usr/local/bin/c3d ' + fnImg + ' -smooth 1mm -resample 50% \-o ' + fnImgRs
        print(cmd)
        os.system(cmd)

    # Preserving data
    poly_data = vtk_read_polydata(fn_mask_ref_vtk)
    print("Mesh data preserved")

    # Propagate Forward
    


    # Propagate Backward



    # Propagate in Full Resolution

"""
read vtk polydata from .vtk file
"""

def vtk_read_polydata(fn):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fn)
    reader.Update()
    return reader.GetOutput()

