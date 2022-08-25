#!/usr/bin/env python3
import argparse
import sys
import json

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
    args = parser.parse_args()

    print(f'img4d: {args.img4d}')
    print(f'seg: {args.seg}')
    print(f'tag: {args.tag}')
    print(f'tpr: {args.tpr}')
    print(f'tpt: {args.tpt}')
    print(f'outdir: {args.outdir}')
    print(f'config: {args.config}')

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
    p.SetUseAffineJitter(True); # remove randomness

    p.Run()


if (__name__ == "__main__"):
    main()
