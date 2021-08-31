import os
from propagation import Propagator

workdir = "/users/jileihao/picsl/propagation/sandbox"

# Create a new Propagator
p = Propagator()

"""bavcta001"""
fnimg = os.path.join(workdir, "test/bavcta001/img4d__bavcta001_trim.nii.gz")
fnseg = os.path.join(workdir, "test/bavcta001/seg03_bavcta001_trim.nii.gz")
fref = 3
targetFrame = [1,2,3,4,5,6,7]

"""bav07"""
# fnimg = os.path.join(workdir, "test/bav07_dcm/bav07.dcm")
# fnseg = os.path.join(workdir, "test/bav07_dcm/seg05_bav07_root_labeled_LPS.nii.gz")
# fref = 5
# targetFrame = [3,5,7]

## Set parameters
p.SetTag("renameTest")
p.SetInputImage(fnimg)
p.SetReferenceSegmentation(fnseg)
p.SetReferenceFrameNumber(fref)
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
p.SetTargetFrames(targetFrame)
p.SetOutputDir(os.path.join(workdir, "out"))
p.SetSmoothingNumberOfIteration(35)
p.SetSmoothingPassband(0.05)

## Add additional mesh to warp
##  parameters: 
##    id: (string) identifier of the mesh. used for list update and deletion, and naming of the file
##    filename: (string) the file path of the mesh
##    smooth: (boolean) indicating if mesh to be smoothed
"""Reference mesh with empty string identifier is added by default and cannot be removed"""
p.AddMeshToWarp('d', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_a.vtk'), True)
p.AddMeshToWarp('e', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_b.vtk'), False)
p.AddMeshToWarp('f', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_c.vtk'), False)

## check list of meshes to be warped
#p.GetWarpingList()
## remove meshe from the list
#p.RemoveMeshFromWarp('c')


### Optional Parameters for testing purpose
#p.SetFullResIterations('5x2x1')
#p.SetDilatedResIteration('5x2x1')
#p.SetGreedyThreads(6)


## Run propagation
"""
    - Set MeshWarpOnly to True to only warp meshes in the WarpingList based on existing 
      registration matrices. 
    - MeshWarpOnly mode is rely on previously generated registration matrices. Therefore, propagation
      with registration has to be completed before MeshWarpOnly run with same:
        - Reference frame
        - Output directory
        - Target frame
        - Tag
      Also do not move, delete, rename any files in the out/tmp folder.
      Or there will be missing file errors.
    - By default MeshWarpOnly is False, meaning propagator will run with full registration
      and mesh warping
"""
p.Run(MeshWarpOnly = True)