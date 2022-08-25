import os
import sys

## An example running the Propagator on windows

sys.path.append('\\Absolute\\path\to\\the\\src');

from src.Propagator import Propagator

workdir = "C:\\playground\\pg_propagation\\Systole\\run_01"

# Create a new Propagator
p = Propagator()

fnimg = "C:\\data\\data_propagation\\Systole\\img4d_seq.nii"
fnseg = "C:\\data\\data_propagation\\Systole\\seg_image.nii.gz"
fref = 1
targetFrame = [2,3,4]


## Set parameters
p.SetTag("winTest")
p.SetInputImage(fnimg)
p.SetReferenceSegmentation(fnseg)
p.SetReferenceFrameNumber(fref)
p.SetGreedyLocation("C:\\tk\greedy\\bin\\greedy.exe")
p.SetVtkLevelSetLocation("C:\\tk\\vtklevelset\\vtklevelset.exe")
p.SetTargetFrames(targetFrame)
p.SetOutputDir(os.path.join(workdir, "out"))
p.SetSmoothingNumberOfIteration(35)
p.SetSmoothingPassband(0.05)

## Add additional mesh to warp
##  parameters: 
##    id: (string) identifier of the mesh. used for list update and deletion, and naming of the file
##    filename: (string) the file path of the mesh
##    smooth: (boolean) indicating if mesh to be smoothed


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
p.Run(MeshWarpOnly = False)