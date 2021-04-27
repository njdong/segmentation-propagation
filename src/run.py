import os
from propagation import Propagator

workdir = "/users/jileihao/playground/sandbox"

# Create a new Propagator
p = Propagator()

# Set Parameters
p.SetTag("dev")
p.SetInputImage(os.path.join(workdir, "test/img4d__bavcta001_trim.nii.gz"))
p.SetReferenceSegmentation(os.path.join(workdir, "test/seg03_bavcta001_trim.nii.gz"))
p.SetReferenceFrameNumber(3)
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
p.SetTargetFrames([1,3,7])
p.SetOutputDir(os.path.join(workdir, "out"))

# Run propagation
p.Run()