import os
import sys

## An example running the propagator on unix-like shell

sys.path.append('/Absolute/path/to/the/src');

from Propagator import Propagator

# Create a new Propagator
p = Propagator()

# Set inputs
onlyWarpMesh = False

## Configure the propagator
p.SetTag("python")
p.SetInputImage("/Users/jileihao/data/data_propagation/greed-2_validation/img4d_seq.nii")
p.SetReferenceSegmentation("/Users/jileihao/data/data_propagation/greed-2_validation/seg_image.nii.gz")
p.SetReferenceFrameNumber(9)
p.SetTargetFrames([10])
p.SetGreedyLocation("/Users/jileihao/dev/greedy-dev/build-release-dynamic/greedy")
p.SetC3dLocation("/Users/jileihao/dev/c3d-dev/build-release-dynamic/c3d")
p.SetVtkLevelSetLocation("/Users/jileihao/dev/cmrep-dev/build-release-dynamic/vtklevelset")
p.SetOutputDir("/Users/jileihao/playground/pg_propagation/greedy-2_validation/python/out")
p.SetSmoothingNumberOfIteration(35)
p.SetSmoothingPassband(0.05)

### Optional Parameters for testing purpose
p.SetFullResIterations('50x50')
p.SetDilatedResIteration('50x50')
p.SetUseAffineJitter(False); # remove randomness
#p.SetGreedyThreads(8)


## Add additional mesh to warp
##  parameters: 
##    id: (string) identifier of the mesh. used for list update and deletion, and naming of the file
##    filename: (string) the file path of the mesh
##    smooth: (boolean) indicating if mesh to be smoothed
"""Reference mesh with empty string identifier is added by default and cannot be removed"""
# p.AddMeshToWarp('d', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_a.vtk'), True)
# p.AddMeshToWarp('e', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_b.vtk'), False)
# p.AddMeshToWarp('f', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_c.vtk'), False)

## check list of meshes to be warped
#p.GetWarpingList()
## remove meshe from the list
#p.RemoveMeshFromWarp('c')




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
p.Run(MeshWarpOnly = onlyWarpMesh)