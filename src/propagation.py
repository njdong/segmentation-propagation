import os
import vtk
import time
import shutil
from LogManager import LogManager
from Image4D import Image4D
from GreedyHelper import GreedyHelper

class Propagator:
    def __init__(self):
        # set default values
        self.fnimg = ""
        self.tag = "default"
        self.fref = -1
        self.outdir = "./out"
        self.greedyLocation = "greedy"
        self.vtklevelsetLocation = "vtklevelset"
        self.greedy = None
        self.multi_res_schedule = '100x100'
        self.metric_spec = 'SSD'

    def SetInputImage(self, _fnimg):
        self.fnimg = _fnimg

    def SetTag(self, _tag):
        self.tag = _tag

    def SetOutputDir(self, _outdir):
        self.outdir = _outdir

    def SetReferenceSegmentation(self, _fnsegref):
        self.fnsegref = _fnsegref
    
    def SetReferenceFrameNumber(self, _fref):
        self.fref = _fref
    
    def SetTargetFrames(self, _fnums):
        """Overrides existing target frames"""
        self.targetFrames = _fnums

    def SetTargetFrameRanges(self, _frange):
        """Overrides existing target frames"""
        self.__parseFrameRangeArray(_frange)

    def SetGreedyLocation(self, _greedyLoc):
        """
            Optional: Set the specific version of greedy for the propagation
            By default the propagator will run greedy from the path
        """
        self.greedyLocation = _greedyLoc

    def SetVtkLevelSetLocation(self, _vtklevelsetLoc):
        """
            Optional: Set the specific version of vtklevelset for the propagation
            By default the propagator will run vtklevelset from the path
        """
        self.vtklevelsetLocation = _vtklevelsetLoc

    def SetFullResIterations(self, _iter):
        """
            Optional: Set the Multi-Resolution Schedule (-n) parameter value for 
            full resolution greedy registration
            Default Value: 100x100
        """
        self.multi_res_schedule = _iter

    def SetMetricSpec(self, _metric_spec):
        """
            Optional: Set the Metric Specification (-m) parameter value for full
            resolution greedy registration
            Default Value: SSD
        """
        self.metric_spec = _metric_spec


    def Run(self):
        self.__propagate()
        print("Propagation completed!")


    def __parseFrameRangeArray(self, rangeArr):
        # todo: replace the placeholder
        self.targetFrames = []



    def __propagate(self):
        #__propagate(fnimg, outdir = "", tag = "", seg_ref = "", framenums = [], fref = 0)

        """
        INPUTS:
        - fnimg: Image filename.
            - Supported formats:
                - (*.dcm) 4D Cartesian DICOM
                - (*.nii) 4D NIfTI
        - outdir: Name of existing directory for output files
        - tag: Study Identifier (e.g., 'bav08_root')
        - seg_ref: Filename of reference segmentation (nii)
        - framenums: List of frame numbers to propagate segementation to
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
        # Initialize Greedy Helper
        self.greedy = GreedyHelper(self.greedyLocation)

        # Validate Input Parameters
        if (self.fnimg == ""):
            raise RuntimeError("Input Image not set!")
        if (len(self.targetFrames) == 0):
            raise RuntimeError("Target Frames not set!")
        if (self.fref == -1):
            raise RuntimeError("Reference Frame Number not set!")
        if (self.fref not in self.targetFrames):
            # always including reference frame in the target frame
            self.targetFrames.append(self.fref)
            self.targetFrames.sort()
            print('Reference frame added to target frames')

        # create output directories
        meshdir = os.path.join(self.outdir, 'mesh')
        if not os.path.exists(self.outdir):
            # directory for output
            os.mkdir(self.outdir)
            # subdirectory for mesh outputs
            os.mkdir(meshdir)

        # recreate mesh dir in case it was removed manually
        if not os.path.exists(meshdir):
            os.mkdir(meshdir)

        tmpdir = os.path.join(self.outdir, "tmp")

        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)

        # performance logger
        perflog = {}
        timepoint = time.time()

        image = None
        
        # parse image type
        if self.fnimg.lower().endswith('.dcm'):
            print('Reading dicom image...')
            # Use the Dicom4D reader to create an Image4D object
            image = Image4D(self.fnimg, 'dicom')
            perflog['Dicom Loading'] = time.time() - timepoint
        elif self.fnimg.lower().endswith(('.nii.gz', '.nii')):
            print('Reading NIfTI image...')
            image = Image4D(self.fnimg, 'nifti')
            # Use the NIfTI reader to create an Image4D object
            perflog['NIfTI Loading'] = time.time() - timepoint
        else:
            print('Unknown image file type')
            return

        

        # Process reference segmentation
        # - Dilate the reference segmentation (mask)
        fn_mask_ref_srs = os.path.join(tmpdir, f'mask_{self.fref}_{self.tag}_srs.nii.gz')
        cmd = f'c3d -int 0 {self.fnsegref} -threshold 1 inf 1 0 \
            -dilate 1 10x10x10vox -resample 50% -o {fn_mask_ref_srs}'
        print("Dilating reference segmentation...")
        print(cmd)
        os.system(cmd)

        # - Create vtk mesh
        fn_mask_ref_vtk = os.path.join(tmpdir, f'mask_{self.fref}_{self.tag}.vtk')
        cmd = f'{self.vtklevelsetLocation} -pl {self.fnsegref} {fn_mask_ref_vtk} 1'
        print("Creating mask mesh...")
        print(cmd)
        os.system(cmd)
        # -- make one copy to the out directory with same naming convention as output mesh
        fn_seg_ref_vtk = os.path.join(meshdir, f'seg_{self.tag}_{self.fref}.vtk')
        shutil.copyfile(fn_mask_ref_vtk, fn_seg_ref_vtk)

        # - Create vtk mesh from the dilated mask
        fn_mask_ref_srs_vtk = os.path.join(tmpdir, f'mask_{self.fref}_{self.tag}_srs.vtk')
        cmd = f"{self.vtklevelsetLocation} -pl {fn_mask_ref_srs} {fn_mask_ref_srs_vtk} 1"
        print("Creating dilated mask mesh...")
        print(cmd)
        os.system(cmd)



        # Export 3D Frames
        # CartesianDicom.Export4D(os.path.join(outdir, 'img4D.nii.gz'))
        # Parallelizable

        framenums = self.targetFrames
        tag = self.tag
        fref = self.fref

        timepoint = time.time()
        
        for i in framenums:
            fnImg = f'{tmpdir}/img_{i}_{tag}.nii.gz'
            fnImgRs = f'{tmpdir}/img_{i}_{tag}_srs.nii.gz'
            image.ExportFrame(i, fnImg)
            
            cmd = 'c3d ' + fnImg + ' -smooth 1mm -resample 50% \-o ' + fnImgRs
            print(cmd)
            os.system(cmd)
        
        perflog['Export 3D Frames'] = time.time() - timepoint
        
        
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

            
            self.__propagation_helper(
                work_dir = tmpdir,
                crnt_ind = i,
                is_forward = True,
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

            self.__propagation_helper(
                work_dir = tmpdir,
                crnt_ind = i,
                is_forward = False,
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
            fn_seg_reslice = f'{self.outdir}/seg_{fref}_to_{fCrnt}_{tag}_reslice.nii.gz'
            fn_seg_reslice_vtk = f'{self.outdir}/mesh/seg_{tag}_{fCrnt}.vtk'

            # transformation filenames
            fn_regout_deform = f'{tmpdir}/warp_{fref}_to_{fCrnt}.nii.gz'
            fn_regout_deform_inv = f'{tmpdir}/warp_{fref}_to_{fCrnt}_inv.nii.gz'

            # run registration and apply warp
            timepoint1 = time.time()

            
            print('Running Full Res Registration...')
            if self.multi_res_schedule != '100x100':
                print(f'Using non-default parameters: -n {self.multi_res_schedule}')
            if self.metric_spec != 'SSD':
                print(f'Using non-default parameter: -m {self.metric_spec}')

            self.greedy.run_reg(
                img_fix = fn_img_fix,
                img_mov = fn_img_mov,
                affine_init = affine_warps,
                regout_deform = fn_regout_deform,
                regout_deform_inv = fn_regout_deform_inv,
                reference_image = fn_reference_frame,
                mask_fix = mask_fix,
                multi_res_schedule = self.multi_res_schedule,
                metric_spec = self.metric_spec
            )
            perflog[f'Full Res Frame {fCrnt} - Registration'] = time.time() - timepoint1
            timepoint1 = time.time()

            print('Applying warp to segmentation...')
            self.greedy.apply_warp(
                image_type = 'label',
                img_fix = fn_img_fix,
                img_mov = self.fnsegref,
                img_reslice = fn_seg_reslice,
                reg_affine = affine_warps,
                reg_deform = fn_regout_deform
            )
            perflog[f'Full Res Frame {fCrnt} - Label Warp'] = time.time() - timepoint1
            timepoint1 = time.time()

            print('Applying warp to mesh...')
            mesh_warps = affine_warps_pts + fn_regout_deform_inv
            self.greedy.apply_warp(
                image_type = 'mesh',
                img_fix = fn_img_fix,
                img_mov = fn_mask_ref_vtk,
                img_reslice = fn_seg_reslice_vtk,
                reg_affine = mesh_warps
            )
            perflog[f'Full Res Frame {fCrnt} - Mesh Warp'] = time.time() - timepoint1
            

        perflog['Full Res Propagation'] = time.time() - timepoint

        fn_perflog = os.path.join(self.outdir, 'perflog.txt')

        if os.path.exists(fn_perflog):
            fLog = open(fn_perflog, 'w')
        else:
            fLog = open(fn_perflog, 'x')

        for k in perflog:
            oneline = f'{k} : {perflog[k]}'
            print(oneline)
            fLog.write(oneline + '\n')

    
    


    def __propagation_helper(self, work_dir, crnt_ind, mask_init, \
        warp_str_array, mask_ref_srs, mask_ref_srs_vtk, is_forward = True):
        """
            Run propagation in one direction (forward or backward), 
            with downsampled images
        """

        fref = self.fref
        framenums = self.targetFrames
        tag = self.tag


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
        self.greedy.run_reg(
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
        self.greedy.apply_warp(
            image_type = 'label',
            img_fix = fn_img_mov,
            img_mov = mask_ref_srs,
            img_reslice = fn_mask_init_reslice,
            reg_affine = warp_str_array[crnt_ind]
        )

        print('Applying warp to mesh...')
        self.greedy.apply_warp(
            image_type = 'mesh',
            img_fix = fn_img_mov,
            img_mov = mask_ref_srs_vtk,
            img_reslice = fn_mask_init_reslice_vtk,
            reg_affine = warp_str_array[crnt_ind]
        )
