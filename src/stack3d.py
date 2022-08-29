#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import SimpleITK as sitk
import os

class Stack3D:
    def __init__(self):
        self._segToStack = {}
        self._tag = ""

    def SetReferenceImage(self, fnImg):
        assert os.path.exists(fnImg), f"Reference image with path '{fnImg}' doesn't exist"
        self._refImg = fnImg
    
    def SetReferenceSegmentation(self, refTP, fnRefSeg):
        assert os.path.exists(fnRefSeg), f"Refernece segmentation image with path '{fnRefSeg}' doesn't exist"
        self._refTP = refTP
        self._segToStack[refTP] = fnRefSeg
    
    def AddSegmentation(self, tp, fnSeg):
        assert os.path.exists(fnSeg), f"Segmentation image with path '{fnSeg}' doesn't exist"
        self._segToStack[tp] = fnSeg

    def SetOutputDir(self, outdir):
        assert os.path.exists(outdir), f"Output directory with path '{outdir}' doesn't exist"
        self._outDir = outdir

    def SetTag(self, tag):
        self._tag = tag

    def Write(self):
        # number of frames
        img_ref = sitk.ReadImage(self._refImg)
        nf = img_ref.GetSize()[3]
        print(f"[Stack3D] Reference image number of frames: {nf}")

        # reference segmentation and blank copy
        seg_ref = sitk.ReadImage(self._segToStack[self._refTP])
        seg_ref_blank = 0*seg_ref

        # append volumes
        vol = []
        print(f"[Stack3D] Start appending volumes")
        for i in range(nf):
            crntTP = i + 1
            print(f"[Stack3D] Processing frame: {crntTP}")
            if crntTP in self._segToStack:
                vol.append(sitk.ReadImage(self._segToStack[crntTP]))
            else:
                vol.append(seg_ref_blank)
            
        # create 4D segmentation series    
        seg4d = sitk.JoinSeries(vol)

        # write 4D segmentation image
        writer = sitk.ImageFileWriter()

        if (self._tag == ""):
            fnout = "seg4d.nii.gz"
        else:
            fnout = f"{self._tag}_seg4d.nii.gz"
            
        writer.SetFileName(os.path.join(self._outDir,fnout))
        writer.Execute(seg4d)