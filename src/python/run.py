import propagation as p
import os

print("CZI-31 python testing script")
workdir = "/users/jileihao/playground/sandbox"
outdir = os.path.join(workdir, "out")
segRefFn = os.path.join(workdir, "seg05_bav07_root_labeled.nii.gz")
dcmFn = os.path.join(workdir, "bav07.dcm")
frameNums = [2, 3, 5, 7, 9, 11]


p.propagate(
    fndcm = dcmFn,
    outdir = outdir,
    tag = "dev",
    seg_ref = segRefFn,
    framenums = frameNums,
    fref = 5
)

print("Completed")