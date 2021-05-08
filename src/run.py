import os
from propagation import Propagator

workdir = "/users/jileihao/playground/sandbox"

# Create a new Propagator
p = Propagator()

# bavcta001
fnimg = os.path.join(workdir, "test/bavcta001/img4d__bavcta001_trim.nii.gz")
fnseg = os.path.join(workdir, "test/bavcta001/seg03_bavcta001_trim.nii.gz")
fref = 3
targetFrame = [1,3]

# bav07
# fnimg = os.path.join(workdir, "test/bav07_dcm/bav07.dcm")
# fnseg = os.path.join(workdir, "test/bav07_dcm/seg05_bav07_root_labeled_LPS.nii.gz")
# fref = 5
# targetFrame = [3,5,7]

# Set Parameters
p.SetTag("renameTest")
p.SetInputImage(fnimg)
p.SetReferenceSegmentation(fnseg)
p.SetReferenceFrameNumber(fref)
p.SetGreedyLocation(os.path.join(workdir, "greedy"))
p.SetVtkLevelSetLocation(os.path.join(workdir, "vtklevelset"))
p.SetTargetFrames(targetFrame)
p.SetOutputDir(os.path.join(workdir, "out"))
### Optional Parameters for testing purpose
p.SetFullResIterations('100x20')
p.SetGreedyThreads(6)


# Run propagation
p.Run()