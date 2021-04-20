import propagation as p
import os

print("CZI-31 python testing script")
workdir = "/users/jileihao/playground/sandbox"
outdir = os.path.join(workdir, "out")
segRefFn = os.path.join(workdir, "seg05_bav07_root_labeled.nii.gz")
fnImg = os.path.join(workdir, "bav07.nii.gz")
frameNums = [1, 2, 3, 4, 5, 6, 7, 8, 9]


p.propagate(
    fnimg = fnImg,
    outdir = outdir,
    tag = "nii",
    seg_ref = segRefFn,
    framenums = frameNums,
    fref = 5
)

print("Completed")