#!/usr/bin/env python

import ipdtools
import argparse
import os
import sys
import logging

def main():

    parser = argparse.ArgumentParser()


    parser.add_argument("--model","-m",
                            help='Choose the model for IPD prediction. See the README of package for more info. DEFAULT: SP2-C2',
                            required=False,
                            choices=["SP2-C2","C2","P4-C2","P5-C3","P6-C4","XL-C2","XL-XL"],
                            default="SP2-C2")

    parser.add_argument("--fastafile","-f",
                            help='Path to a fasta file.',
                            required=True,
                            default=None)

    parser.add_argument('--output_csv',"-o",
                            help='Output CSV file of predicted IPDs.',
                            required=True)

    parser.add_argument('--verbosity',"-v",
                            help='Choose your verbosity. Default: INFO',
                            required=False,
                            default="WARNING",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR, CRITICAL"])

    args = parser.parse_args()

    fastafile = os.path.realpath(args.fastafile)
    output_csv = os.path.realpath(args.output_csv)

    verboselevel = "logging."+str(args.verbosity)
    logging.basicConfig(stream=sys.stdout, level=eval(verboselevel),
                        format='%(asctime)s %(message)s')

    show_progress_bar = False
    if args.verbosity in ["DEBUG","INFO"]:
        show_progress_bar = True

    ipdtools.ipdModel.compute_fasta_to_csv(modelname=args.model,fastafile=fastafile,csvout=output_csv,show_progress_bar=show_progress_bar)

if __name__ == "__main__":
    main()
