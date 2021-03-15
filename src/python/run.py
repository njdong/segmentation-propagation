import propagation as p
import os

print("CZI-31 python testing script")
workdir = "/users/jileihao/playground/sandbox"
outdir = os.path.join(workdir, "out")
segRefFn = os.path.join(workdir, "seg05_bav07_root_labeled.nii.gz")
dcmFn = os.path.join(workdir, "bav07.dcm")
frameNums = [1, 5, 6, 7]


p.propagate( \
    fndcm = dcmFn, \
    outdir = outdir, \
    tag = "dev", \
    seg_ref = segRefFn, \
    framenums = frameNums, \
    fref = 5 \
)


"""
import vtk_read_polydata as vrp
meshFn = os.path.join(workdir, "seg05_bav07_root_labeled.vtk")
reader = vrp.vtk_read_polydata(meshFn)
poly_data = reader.GetOutput()
#print(poly_data)
print(poly_data.GetNumberOfCells())
"""

"""
print("POLY DATA --------------------------------------")
print(poly_data)
print("STRIPS -------------------------------------")
strips = poly_data.GetStrips()
print(strips)

# Restore cell data back to poly data
poly_data.SetStrips(strips)
print("RESULT -----------------------------------------")
print(poly_data)

writer = vtk.vtkPolyDataWriter()
writer.SetFileName(os.path.join(workdir, "mesh_out.vtk"))
writer.SetInputData(poly_data)
writer.Write()
"""




print("Complete")