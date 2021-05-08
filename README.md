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


p.SetTargetFrames(targetFrame)
p.SetOutputDir(os.path.join(workdir, "out"))

# Optional Parameters for testing purpose
# -- Set to use specific version of greedy.
# -- If not set, it will use greedy in system PATH
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
# -- Set to use specific version of vtklevelset.
# -- If not set, it will use greedy in system PATH
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
# -- Set the iteration schedule for the full resolution propagation
# -- If not set, it will use 100x100
p.SetFullResIterations('100x20')
# -- Set number of threads for greedy to use
# -- If not set, it will use all threads available
p.SetGreedyThreads(6)

# Run propagation
p.Run()
```


## Version History
### Version 1.0
- czi-33: First usable version
### Version 1.1
- czi-67: Modfied to adopt latest greedy and vtkleveset
### Version 1.2
- czi-73: Reorganizes mesh output files into a dedicated mesh folder with new naming convention
- czi-74: Added mesh point data renaming logic; Added multiple options for configure greedy.
