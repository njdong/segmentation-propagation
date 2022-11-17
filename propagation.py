#!/usr/bin/env python3
import json
import sys
import argparse

sys.path.append("./src")

from src.Propagator import Propagator

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("img4d", help="Filepath to the 4D input image")
    parser.add_argument("seg", help="Filepath to the reference segmentation image")
    parser.add_argument("tag", help="Identifier for current run")
    parser.add_argument("tpr", help="Reference timepoint", type=int)
    parser.add_argument("tpt", help="Semicolon separated string list of target timepoints")
    parser.add_argument("outdir", help="Output directory for temp files and results")
    parser.add_argument("config", help="Filepath to a config json file")

    parser.add_argument("-warp_only", action='store_true',
                        help="If specified, run warp without registration")
    parser.add_argument("-add_mesh", nargs=2, action="append",
                        help="Add an additional mesh to warp. Format: -add_mesh tag filename")
    
    args = parser.parse_args()

    print(f'img4d: {args.img4d}')
    print(f'seg: {args.seg}')
    print(f'tag: {args.tag}')
    print(f'tpr: {args.tpr}')
    print(f'tpt: {args.tpt}')
    print(f'outdir: {args.outdir}')
    print(f'config: {args.config}')
    print(f'-warp_only: {args.warp_only}')
    print(f'-add_mesh: {args.add_mesh}')

    p = Propagator()

    # Configure the propagator with args
    p.SetTag(args.tag)
    p.SetInputImage(args.img4d)
    p.SetReferenceSegmentation(args.seg)
    p.SetReferenceFrameNumber(args.tpr)
    p.SetOutputDir(args.outdir)

    tptParsed = list(map(int, args.tpt.split(';')))
    print(f'parsed tpt= {tptParsed}')
    p.SetTargetFrames(tptParsed)

    f = open(args.config)
    config = json.load(f)

    # Read from configuration file
    p.SetGreedyLocation(config["greedyLocation"])
    p.SetC3dLocation(config["c3dLocation"])
    p.SetVtkLevelSetLocation(config["vtklevelsetLocation"])
    p.SetFullResIterations(config["fullResIteration"])
    p.SetDilatedResIteration(config["downSampledIteration"])
    p.SetSmoothingNumberOfIteration(int(config["meshSmoothingNumberOfIterations"]))
    p.SetSmoothingPassband(float(config["meshSmoothingPassband"]))
    p.SetDisablePythonMesh(config["disablePythonMesh"])
    p.SetUseAffineJitter(config["useAffineJitter"])  # remove randomness

    if args.add_mesh:
        for mesh_input in args.add_mesh:
            p.AddMeshToWarp(mesh_input[0], mesh_input[1], True)

    p.Run(warpOnly=args.warp_only)


if (__name__ == "__main__"):
    main()
