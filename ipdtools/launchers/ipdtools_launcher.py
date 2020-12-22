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
                            help='Choose the model for IPD prediction. See the package\'s README for more info. DEFAULT: SP2-C2 (PacBio Sequel I)',
                            required=False,
                            choices=["SP2-C2","C2","P4-C2","P5-C3","P6-C4","XL-C2","XL-XL"],
                            default="SP2-C2")

    parser.add_argument("--fastafile","-f",
                            help='Path to a fasta file.',
                            required=True,
                            default=None)

    parser.add_argument('--output_csv',"-o",
                            help='Output CSV file of predicted IPDs. See the README for further details on the output.',
                            required=True)

    parser.add_argument('--verbosity',"-v",
                            help='Choose your verbosity. Default: INFO',
                            required=False,
                            default="INFO",
                            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    parser.add_argument('--progress_bar',"-p",
                            help='Displays a progress bar',
                            action='store_true')

    parser.add_argument('--nproc','-n',
                        help="Max number of processors for parallelism. DEFAULT: 1, positive integers only. Rq: The programm will actually use n+1 >= 2 CPU no matter what.",
                        required=False,
                        default=1,
                        type=check_positive)

    parser.add_argument('--indexing','-i',
                        help="Is the indexing corresponding to the .fasta reference 1-based (default PacBio) or 0-based. DEFAULT: 1",
                        required=False,
                        default=1,
                        type=int,
                        choices=[0,1])

    args = parser.parse_args()

    fastafile = os.path.realpath(args.fastafile)
    output_csv = os.path.realpath(args.output_csv)

    verboselevel = "logging."+str(args.verbosity)
    logging.basicConfig(level=eval(verboselevel),
                        format='%(asctime)s %(message)s',
                        stream=sys.stdout)


    show_progress_bar = False
    if args.progress_bar:
        show_progress_bar = True

    ipdtools.ipdModel.compute_fasta_to_csv(modelname=args.model, fastafile=fastafile, csv_out=output_csv,
                                           show_progress_bar=show_progress_bar,nproc=args.nproc,indexing=args.indexing)

if __name__ == "__main__":
    main()
