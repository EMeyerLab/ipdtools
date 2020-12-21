#!/usr/bin/env python

import ipdtools
import argparse
import os
import sys
import logging

def main():

    parser = argparse.ArgumentParser()

    def check_positive(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue

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

    parser.add_argument('--nproc','-n',
                        help="Max number of processors for parallelism. DEFAULT: 1 [NOT IMPLEMENTED YET]",
                        required=False,
                        default=1,
                        type=check_positive)


    args = parser.parse_args()

    fastafile = os.path.realpath(args.fastafile)
    output_csv = os.path.realpath(args.output_csv)

    verboselevel = "logging."+str(args.verbosity)
    logging.basicConfig(stream=sys.stdout, level=eval(verboselevel),
                        format='%(asctime)s %(message)s')

    show_progress_bar = False
    if args.verbosity in ["DEBUG","INFO"]:
        show_progress_bar = True

    ipdtools.ipdModel.compute_fasta_to_csv(modelname=args.model,fastafile=fastafile,csvout=output_csv,show_progress_bar=show_progress_bar,nproc=args.nproc)

if __name__ == "__main__":
    main()
