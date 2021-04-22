import os
from propagation import Propagator

workdir = "/users/jileihao/playground/sandbox"

# Create a new Propagator
p = Propagator()

# Set Parameters
p.SetTag("dev")
p.SetInputImage(os.path.join(workdir, "bav07.nii.gz"))
p.SetReferenceSegmentation(os.path.join(workdir, "seg05_bav07_root_labeled.nii.gz"))
p.SetReferenceFrameNumber(5)
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
p.SetTargetFrames([1,5,8])
p.SetOutputDir(os.path.join(workdir, "out"))

# Run propagation
p.Run()