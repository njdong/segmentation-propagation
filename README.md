# segmentation-propagation

The function propagation.m is used to transform a segmentation in one frame to other frames within an image series. Deformable registration between frames is performed with the "greedy" tool.

## Example
```python
import os
from propagation import Propagator


# Optional: Recommend a work directory to organize all inputs and outputs
workdir = "/users/guest/playground"

# Create a new Propagator
p = Propagator()


fnimg = os.path.join(workdir, "img4d.nii.gz")
fnseg = os.path.join(workdir, "seg03.nii.gz")
fref = 3
targetFrame = [1,3,7]

# Configure Propagation Parameters
# -- Set an identifier for current run.
# -- Do not include illegal characters for a file name
p.SetTag("Test")

# -- Input Image filename
p.SetInputImage(fnimg)

# -- Reference segmentation is the segmentation to be propagated to target frames
p.SetReferenceSegmentation(fnseg)

# -- Reference frame number is the frame number of the segmentation to be propagated
p.SetReferenceFrameNumber(fref)

# -- Target frames for propagation
p.SetTargetFrames(targetFrame)

# -- Output directory for results and itermediate files
p.SetOutputDir(os.path.join(workdir, "out"))

## Add additional mesh to warp
##  parameters: 
##    id: (string) identifier of the mesh. used for list update and deletion, and naming of the file
##    filename: (string) the file path of the mesh
##    smooth: (boolean) indicating if mesh to be smoothed
"""Reference mesh with empty string identifier is added by default and cannot be removed"""
p.AddMeshToWarp('a', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_a.vtk'), True)
p.AddMeshToWarp('b', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_b.vtk'), False)
p.AddMeshToWarp('c', os.path.join(workdir, 'test/bavcta001/seg03_bavcta001_c.vtk'), False)

## check list of meshes to be warped
#p.GetWarpingList()
## remove meshe from the list
#p.RemoveMeshFromWarp('c')

### Optional Parameters for debugging purpose
#p.SetFullResIterations('5x2x1')
#p.SetDilatedResIteration('5x2x1')
#p.SetGreedyThreads(6)

# Optional Parameters for testing purpose
# -- Set to use specific version of greedy.
# -- By default, it will use greedy in the system PATH
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
# -- Set to use specific version of vtklevelset.
# -- By default, it will use greedy in the system PATH
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
# -- Set the iteration schedule for the full resolution propagation
# -- By default, it will use 100x100
p.SetFullResIterations('100x20')
# -- Set number of threads for greedy to use
# -- By default, it will use all threads available
p.SetGreedyThreads(6)

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
```


## Version History
### Version 1.0
- czi-33: First usable version
### Version 1.1
- czi-67: Modfied to adopt latest greedy and vtkleveset
### Version 1.2
- czi-73: Reorganizes mesh output files into a dedicated mesh folder with new naming convention
- czi-74: Added mesh point data renaming logic; Added multiple options for configure greedy.
- czi-76: Added Mesh Warping feature; Added run mode that can warp multiple meshes using existing
          registration matrices without running registration again.
### Version 1.2.1
- czi-105: Added flag to control smoothing for added mesh to warp