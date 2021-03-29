import os
import vtk
import time
import asyncio
from Dicom4D import Dicom4D



def propagate(fnimg, outdir = "", tag = "", seg_ref = "", framenums = [], fref = 0):

    """
    INPUTS:
    - fnimg: Image filename.
        - Supported formats:
            - (*.dcm) 4D Cartesian DICOM
            - (*.nii) 4D NIfTI
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

    if len(framenums) == 0:
        return

    # directory storing temporary files
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tmpdir = os.path.join(outdir, "tmp")

    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    # performance logger
    perflog = {}
    timepoint = time.time()
    
    # parse image type
    if fnimg..lower().endswith('.dcm'):
        print('Reading dicom image...')
        # Use the Dicom4D reader to create an Image4D object
        CartesianDicom = Dicom4D(fndcm)
        perflog['Dicom Loading'] =  time.time() - timepoint
    elif fnimg.lower().endswith(('.nii.gz', '.nii')):
        print('Reading NIfTI image...')
        # Use the NIfTI reader to create an Image4D object


    # Process reference segmentation
    # - Dilate the reference segmentation (mask)
    fn_mask_ref_srs = os.path.join(tmpdir, f'mask_{fref}_{tag}_srs.nii.gz')
    cmd = f'c3d -int 0 {seg_ref} -threshold 1 inf 1 0 \
        -dilate 1 10x10x10vox -resample 50% -o {fn_mask_ref_srs}'
    print("Dilating reference segmentation...")
    print(cmd)
    os.system(cmd)

    # - Create vtk mesh
    fn_mask_ref_vtk = os.path.join(tmpdir, f'mask_{fref}_{tag}.vtk')
    cmd = f'vtklevelset {seg_ref} {fn_mask_ref_vtk} 1'
    print("Creating mask mesh...")
    print(cmd)
    os.system(cmd)

    # - Create vtk mesh from the dilated mask
    fn_mask_ref_srs_vtk = os.path.join(tmpdir, f'mask_{fref}_{tag}_srs.vtk')
    cmd = f"vtklevelset {fn_mask_ref_srs} {fn_mask_ref_srs_vtk} 1"
    print("Creating dilated mask mesh...")
    print(cmd)
    os.system(cmd)



    # Export 3D Frames
    # CartesianDicom.Export4D(os.path.join(outdir, 'img4D.nii.gz'))
    # Parallelizable

    timepoint = time.time()
    
    for i in framenums:
        fnImg = f'{tmpdir}/img_{i}_{tag}.nii.gz'
        fnImgRs = f'{tmpdir}/img_{i}_{tag}_srs.nii.gz'
        CartesianDicom.ExportFrame(i, fnImg)
        
        cmd = 'c3d ' + fnImg + ' -smooth 1mm -resample 50% \-o ' + fnImgRs
        print(cmd)
        os.system(cmd)
    
    perflog['Export 3D Frames'] = time.time() - timepoint

    
    # Preserving data
    polyData = vtk_read_polydata(fn_mask_ref_vtk)
    print("Mesh data preserved")
    

    # Initialize warp string array
    warp_str_array = [''] * len(framenums)

    ref_ind = framenums.index(fref)

    
    # Propagate Forward
    print('---------------------------------')
    print('Propagating forward')
    print('---------------------------------')

    timepoint = time.time()

    for i in range(ref_ind, len(framenums) - 1):
        fCrnt = framenums[i]
        fNext = framenums[i + 1]
        fPrev = framenums[i - 1]

        print('Current Frame: ', framenums[i])

        if fCrnt == fref:
            fn_mask_init = fn_mask_ref_srs
        else:
            fn_mask_init = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'

        
        propagation_helper(
            work_dir = tmpdir,
            framenums = framenums,
            tag = tag,
            crnt_ind = i,
            is_forward = True,
            fref = fref,
            mask_init = fn_mask_init,
            warp_str_array = warp_str_array,
            mask_ref_srs = fn_mask_ref_srs,
            mask_ref_srs_vtk = fn_mask_ref_srs_vtk)
        
    perflog['Forward Propagation'] = time.time() - timepoint
    timepoint = time.time()

    # Propagate Backward
    print('---------------------------------')
    print('Propagating backward')
    print('---------------------------------')

    # - Clean up warp str array
    warp_str_array = [''] * len(framenums)

    for i in range(ref_ind, 0, -1):
        fCrnt = framenums[i]
        fNext = framenums[i - 1]
        fPrev = framenums[i + 1] if fCrnt != fref else -1

        print('Current Frame: ', framenums[i])

        if fCrnt == fref:
            fn_mask_init = fn_mask_ref_srs
        else:
            fn_mask_init = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'

        propagation_helper(
            work_dir = tmpdir,
            framenums = framenums,
            tag = tag,
            crnt_ind = i,
            is_forward = False,
            fref = fref,
            mask_init = fn_mask_init,
            warp_str_array = warp_str_array,
            mask_ref_srs = fn_mask_ref_srs,
            mask_ref_srs_vtk = fn_mask_ref_srs_vtk)
    
    perflog['Backward Propagation'] = time.time() - timepoint
    timepoint = time.time()

    # Propagate in Full Resolution
    
    print('---------------------------------')
    print('Propagating in Full Resolution: ')
    print('---------------------------------')

    for i in range(0, len(framenums)):
        fCrnt = framenums[i]
        print('Processing frame: ', fCrnt)
        
        if fCrnt == fref:
            continue

        fn_img_fix = f'{tmpdir}/img_{fCrnt}_{tag}.nii.gz'
        fn_img_mov = f'{tmpdir}/img_{fref}_{tag}.nii.gz'

        # recursively build affine warp
        affine_warps = ''
        affine_warps_pts = ''

        if i > ref_ind:
            fPrev = framenums[i - 1]

            for j in range (ref_ind + 1, i + 1):
                fn_regout_affine_init = f'{tmpdir}/affine_{framenums[j]}_to_{framenums[j - 1]}_srs_init.mat'
                affine_warps = affine_warps + ' ' + fn_regout_affine_init + ',-1 '
                affine_warps_pts = fn_regout_affine_init + ' ' + affine_warps_pts
        else:
            fPrev = framenums[i + 1]

            for j in range (ref_ind - 1, i - 1, -1):
                fn_regout_affine_init = f'{tmpdir}/affine_{framenums[j]}_to_{framenums[j + 1]}_srs_init.mat'
                affine_warps = affine_warps + ' ' + fn_regout_affine_init + ',-1 '
                affine_warps_pts = fn_regout_affine_init + ' ' + affine_warps_pts

        print("affine_warps: ", affine_warps)
        print("affine_warps_pts: ", affine_warps_pts)

        
        # full resolution mask for this frame
        mask_fix_srs = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_srs_reslice_init.nii.gz'
        mask_fix = f'{tmpdir}/mask_{fPrev}_to_{fCrnt}_{tag}_reslice_init.nii.gz'
        print('Generating full res mask...')
        cmd = f'c3d -interpolation NearestNeighbor {fn_img_fix} {mask_fix_srs} -reslice-identity -o {mask_fix}'
        print(cmd)
        os.system(cmd)

        # trim to generate a reference frame
        fn_reference_frame = f'{tmpdir}/reference_{fPrev}_to_{fCrnt}_{tag}.nii.gz'
        cmd = f'c3d {mask_fix} -trim 0vox -o {fn_reference_frame}'
        print(cmd)
        os.system(cmd)

        
        # output file location
        fn_seg_reslice = f'{outdir}/seg_{fref}_to_{fCrnt}_{tag}_reslice.nii.gz'
        fn_seg_reslice_vtk = f'{outdir}/seg_{fref}_to_{fCrnt}_{tag}_reslice_labeled.vtk'

        # transformation filenames
        fn_regout_deform = f'{tmpdir}/warp_{fref}_to_{fCrnt}.nii.gz'
        fn_regout_deform_inv = f'{tmpdir}/warp_{fref}_to_{fCrnt}_inv.nii.gz'

        # run registration and apply warp
        timepoint1 = time.time()

        
        print('Running Full Res Registration...')
        greedy_call(
            img_fix = fn_img_fix,
            img_mov = fn_img_mov,
            affine_init = affine_warps,
            regout_deform = fn_regout_deform,
            regout_deform_inv = fn_regout_deform_inv,
            reference_image = fn_reference_frame,
            mask_fix = mask_fix
        )
        perflog[f'Full Res Frame {fCrnt} - Registration'] = time.time() - timepoint1
        timepoint1 = time.time()

        print('Applying warp to segmentation...')
        apply_warp(
            image_type = 'label',
            img_fix = fn_img_fix,
            img_mov = seg_ref,
            img_reslice = fn_seg_reslice,
            reg_affine = affine_warps,
            reg_deform = fn_regout_deform
        )
        perflog[f'Full Res Frame {fCrnt} - Label Warp'] = time.time() - timepoint1
        timepoint1 = time.time()

        print('Applying warp to mesh...')
        mesh_warps = affine_warps_pts + fn_regout_deform_inv
        apply_warp(
            image_type = 'mesh',
            img_fix = fn_img_fix,
            img_mov = fn_mask_ref_vtk,
            img_reslice = fn_seg_reslice_vtk,
            reg_affine = mesh_warps
        )
        perflog[f'Full Res Frame {fCrnt} - Mesh Warp'] = time.time() - timepoint1

        # restore cell data 
        print('Restoring poly data to mesh...')
        vtk_replace_stripes(fn_seg_reslice_vtk, polyData)
        

    perflog['Full Res Propagation'] = time.time() - timepoint

    fn_perflog = os.path.join(outdir, 'perflog.txt')

    if os.path.exists(fn_perflog):
        fLog = open(fn_perflog, 'w')
    else:
        fLog = open(fn_perflog, 'x')

    for k in perflog:
        oneline = f'{k} : {perflog[k]}'
        print(oneline)
        fLog.write(oneline + '\n')

    


def propagation_helper(work_dir, framenums, tag, crnt_ind, fref, mask_init, \
    warp_str_array, mask_ref_srs, mask_ref_srs_vtk, is_forward = True):
    # forward propagation is in incremental order, backward is the reverse
    next_ind = crnt_ind + 1 if is_forward else crnt_ind - 1

    fCrnt = framenums[crnt_ind]
    fNext = framenums[next_ind]
    
    # current dilated image as fix
    fn_img_fix = f'{work_dir}/img_{fCrnt}_{tag}_srs.nii.gz'
    # next dilated image as moving image
    fn_img_mov = f'{work_dir}/img_{fNext}_{tag}_srs.nii.gz'
    
    # filenames of initial transformations
    fn_regout_affine = f'{work_dir}/affine_{fNext}_to_{fCrnt}_srs_init.mat'
    fn_regout_deform = f'{work_dir}/warp_{fNext}_to_{fCrnt}_srs_init.nii.gz'
    fn_regout_deform_inv = f'{work_dir}/warp_{fNext}_to_{fCrnt}_srs_init_inv.nii.gz'
    
    # call greedy to generate transformations
    greedy_call(
        img_fix = fn_img_fix,
        img_mov = fn_img_mov,
        regout_affine = fn_regout_affine,
        regout_deform = fn_regout_deform,
        regout_deform_inv = fn_regout_deform_inv,
        mask_fix = mask_init
    )

    # Build warp string array recursively
    if fCrnt == fref:
        warp_str_array[crnt_ind] = f'{fn_regout_affine},-1 {fn_regout_deform_inv} '
    else:
        prev_ind = crnt_ind - 1 if is_forward else crnt_ind + 1
        warp_str_array[crnt_ind] = f'{warp_str_array[prev_ind]} {fn_regout_affine},-1 {fn_regout_deform_inv} '

    # Parallelizable
    fn_mask_init_reslice = f'{work_dir}/mask_{fCrnt}_to_{fNext}_{tag}_srs_reslice_init.nii.gz'
    fn_mask_init_reslice_vtk = f'{work_dir}/mask_{fCrnt}_to_{fNext}_{tag}_srs_reslice_init.vtk'

    # call greedy applying warp
    print('Applying warp to label...')
    apply_warp(
        image_type = 'label',
        img_fix = fn_img_mov,
        img_mov = mask_ref_srs,
        img_reslice = fn_mask_init_reslice,
        reg_affine = warp_str_array[crnt_ind]
    )

    print('Applying warp to mesh...')
    apply_warp(
        image_type = 'mesh',
        img_fix = fn_img_mov,
        img_mov = mask_ref_srs_vtk,
        img_reslice = fn_mask_init_reslice_vtk,
        reg_affine = warp_str_array[crnt_ind]
    )

def vtk_replace_stripes(fnmesh, polyData):
    """
    Replace mesh with polyData stripes
    """
    data = vtk_read_polydata(fnmesh)
    data.SetStrips(polyData.GetStrips())
    vtk_write_polydata(fnmesh, data)

def vtk_write_polydata(fn, data):
    """
    Write data to vtk mesh file
    """
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fn)
    writer.SetInputData(data)
    writer.Write()

def vtk_read_polydata(fn):
    """
    Read vtk polydata from .vtk file
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(fn)
    reader.Update()
    return reader.GetOutput()


def greedy_call(img_fix, img_mov, regout_deform_inv, mask_fix, \
    affine_init = '', regout_affine = '', regout_deform = '', reference_image = ''):

    """
    Make system call to greedy

    Calls greedy for deformable registration between two images. Creates
    affine transformation file and/or deformation field depending on which
    optional input filenames are provided.

    INPUT:
    - img_fix: fixed (reference) image filename
    - img_mov: moving image filename
    - affine_init: filename of affine transform for initialization (optional)
    - regout_affine: filename of output affine transformation (optional)
    - regout_deform: filename of output deformation (optional)
    - regout_deform_inv: filename of output inverse deformation (optional)
    - mask_fix: filename of mask for the fixed image
    
    The optional input arguments determine which parameters are used in the 
    call to greedy (affine, deformable registration, or both). 

    Use full paths for all filenames.
    """

    cmd = ''

    if regout_affine != '' and regout_deform != '':
        # Affine generation
        cmd = f'greedy -d 3 \
            -a \
            -i {img_fix} {img_mov} '

        if reference_image != '':
            cmd = cmd + f'-rf {reference_image} '
        
        cmd = cmd + f'-ia-identity \
            -dof 6 \
            -s 3mm 1.5mm \
            -gm {mask_fix} \
            -o {regout_affine} '
        
        print('greedy_call: ', cmd)
        os.system(cmd)

        # Deform generation
        cmd = f'greedy -d 3 \
            -i {img_fix} {img_mov} \
            -it {regout_affine} '

        if reference_image != '':
            cmd = cmd + f'-rf {reference_image} '

        cmd = cmd + f'-m SSD \
            -n 100x100 \
            -s 3mm 1.5mm \
            -gm {mask_fix} \
            -o {regout_deform} '

        if regout_deform_inv != '':
            cmd = cmd + f' -oinv {regout_deform_inv} '

        print('greedy_call: ', cmd)
        os.system(cmd)
        print('greedy_call: Affine + Deformable registration computed!')
    elif regout_affine == '' and regout_deform != '':
        if regout_deform_inv != '':
            if affine_init != '':
                cmd = f'greedy -d 3 \
                    -i {img_fix} {img_mov} '
                if reference_image != '':
                    cmd = cmd + f'-rf {reference_image} '
                
                cmd = cmd + f'-m SSD \
                    -n 100x100 \
                    -it {affine_init} \
                    -gm {mask_fix} \
                    -s 3mm 1.5mm \
                    -o {regout_deform} \
                    -oinv {regout_deform_inv}'
                print('Initialized with input affine transform')
                print('greedy_call: ', cmd)
                os.system(cmd)
                print('greedy_call: Only deformable registration computed!')
        
                
    
    



def apply_warp(image_type, img_fix, img_mov, img_reslice, \
    reg_affine = '', reg_deform = '', reg_deform_inv = ''):
    """
    Calls greedy to apply an affine and/or other deformation to a grayscale
    image, label map (segmentation), or vtk mesh. 

    INPUT:
    - image_type: 'grayscale', 'label', or 'mesh'
    - img_fix: fixed (reference) image filename
    - img_reslice: filename of output (warped) image
    - reg_affine: filename of affine registration (optional)
    - reg_deform: filename of deformation field (optional)
    - reg_deform_inv: filename of inverse deformation field (optional)
    
    The optional input arguments determine which parameters are used in the 
    call to greedy (affine, deformable registration, or both). 
    
    Use full paths for all filenames.
    """
    
    cmd = f'greedy -d 3 \
        -rf {img_fix} \
        -r {reg_deform} {reg_affine} '
    
    if image_type == 'grayscale':
        cmd = cmd + f' \
            -ri LINEAR \
            -rm {img_mov} {img_reslice} '
    elif image_type == 'label':
        cmd = cmd + f' \
            -ri NN \
            -rm {img_mov} {img_reslice} '
    elif image_type == 'mesh':
        cmd = cmd + f' \
            -rs {img_mov} {img_reslice} '

    print('apply_warp: ', cmd)
    os.system(cmd)
    print('apply_warp: Affine + Deformable transformation applied!')

