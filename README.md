The PyPi release is still in alpha.
Automated testings, docs and other things are still missing.
Please contact me if you can't use the package using the following instructions

# ipdtools

Retro-ingineering of PacBio's modelling of IPD (**RS II, Sequel I, Sequel II v2**). This **won't** work for **Sequel II with v1.0 chemistry.**

# Background

The PacBio's SMRT sequencing allows the measurement of the DNA polymerase's kinetics. It has been shown that slowing of the polymerase can be linked to:

- DNA Modifications
- Secondary structures

This kinetic can then be compared with the kinetic of a control sample from a **Whole Genome Amplification** (WGA) to detect possible modifications of nucleotides within DNA:

- **N6-methyl-adenine (n6mA)**
- **5-methyl-cytosine (5mC)**
- **4-methyl-cytosine (4mC)**

As it is known that the polymerase's speed at a **nucleotide N** is influenced by the surrounding **N-5/N+4 nucleotides**, it is most of the times possible to **predict** the polymerase's speed if the given nucleotide is methylated or not.

  For people who don't use WGA, PacBio provides a **modelling of the kinetics of the polymerase** for unmethylated DNA and methylated DNA, that uses **Machine-Learning** of the **surrounding context** of nucleotides to predict the **Inter-Pulse Duration** (IPD), which is directly the speed of the polymerase.

  **Unfortunately, this ML modeling is not directly accessible to scientists**

> From the source code of **kineticsTools** and the **models** published by PacBio at [kineticsTools](https://github.com/PacificBiosciences/kineticsTools), it was possible to make a **user-friendly** command-line-interface (**CLI**) and **python API** to easily predict **any IPD, in any context** for unmethylated DNA.

NB: Methylated DNA is harder to predict because it often involves a mixture of many possible combinations of signals.


# Requirements

- Python >= 3.6

- Any C/C++ compiler correctly installed

> Works well on Ubuntu 18.04 LTS x86. 
> 
> Hasn't been tested on either Windows or macOS (especially not the last ARM-based versions of MacOS)

# Installation

ipdtools can be installed with **pip**, which is automatically distributed with any version of python

```bash
git clone https://github.com/GDelevoye/ipdtools.git
pip install ./ipdtools
```

Please note that, due to the fact we're compiling a third-party dependency, the "editable" pip install will not be available. Do not use pip install -e ./ipdtools

# Usage

```console
(base) user@computer:~/ipdtools/resources$ ipdtools -f test.fasta -o test.csv -p -n 6 --model P6-C4 --verbosity INFO
2020-12-22 15:14:29,067 [INFO] Using 6 CPUs
2020-12-22 15:14:29,079 [INFO] Number of nucleotides to handle : 463
Nucleotides: 100%|███████████████████████████████| 463/463 [00:10<00:00, 39.59it/s]
2020-12-22 15:14:39,145 [INFO] Result saved at /export/home1/users/gmc/delevoye/Github/ipdtools/ipdtools/resources/test.csv
2020-12-22 15:14:39,148 [INFO] DONE
```

## Command-line interface (CLI)

### Generate a .csv of IPDs from a fasta file:

```
user@computer$:ipdtools -f myfile.fasta -o output.csv -m [MODELNAME]
```

> **By default, [MODELNAME]  is SP2-C2**, which is supposed to work on any **Sequel I** or **Sequel IIv2** Data. You might want to change it if you work with older **RSII** Data.

> See ["Which model should I use ?"](#whichmodel) for more info

Detailed Usage

```console
user@computer$:ipdtools --help
usage: ipdtools [-h] [--model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL}]
                --fastafile FASTAFILE --output_csv OUTPUT_CSV
                [--verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                [--progress_bar] [--nproc NPROC] [--indexing {0,1}]

optional arguments:
  -h, --help            show this help message and exit
  --model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL}, -m {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL}
                        Choose the model for IPD prediction. See the README of
                        package for more info. DEFAULT: SP2-C2
  --fastafile FASTAFILE, -f FASTAFILE
                        Path to a fasta file.
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Output CSV file of predicted IPDs.
  --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Choose your verbosity. Default: INFO
  --progress_bar, -p    Displays a progress bar
  --nproc NPROC, -n NPROC
                        Max number of processors for parallelism. DEFAULT: 1,
                        positive integers only. Rq: The programm will actually
                        use n+1 >= 2 CPU no matter what.
  --indexing {0,1}, -i {0,1}
                        Is the indexing corresponding to the .fasta reference
                        1-based (default PacBio) or 0-based
```

### Output .csv header

The output will be a .csv file like this one:

Fasta_ID | Position | Strand | Nucleotide | Prediction
:-: |:-: | :-: | :-: | :-: |
[CONTIG] | 0 | `0` | A | xxxx
[CONTIG] | 0 | `1` | T | xxxx
[CONTIG] | 1 | `0` | T | xxxx
[CONTIG] | 1 | `1` | A | xxxx
[CONTIG] | 2 | `0` | C | xxxx
[CONTIG] | 2 | `1` | G | xxxx
[CONTIG] | 3 | `0` | G | xxxx
[CONTIG] | 3 | `1` | C | xxxx
[CONTIG] | 4 | `0` | A | xxxx
[CONTIG] | 4 | `1` | T | xxxx


Of course, depending on using a 0-based or 1-based indexing, the columns will vary.

> The csv is sorted (in ascending order) by:
> * Fasta_ID
> * Position
> * Strand

## Python API

ipdtools can be used in python to predict from a fasta file


it can also be used to predict directly from a python string with **Str2IPD**:

```python
>>> import ipdtools
>>> predictor = ipdtools.ipdModel.Str2IPD("ATGCTAGCTTTTTGNCTGATTAGCTGA",model="SP2-C2",indexing=0) # Default model is SP2-C2
>>> predictor.predict(position=0) # The default is strand0
1.0641327
>>> predictor.predict(position=0,strand=1)
0.42285475
```

**WARNING** : This feature is experimental. PacBio's model can produce outputs even for very incoherent parameters (like, "Strand 4, Position -1") without raising any error or warning. Although this behaviour is really not reasuring at first, I made sure with manual checkings that the predictions made from valid parameters stay valid (they apparently do). So, when using correct parameters, there is no problem. The problem is only there IF you're using wrong parameters (like, seeking the ipds at position 350 if your sequence is 349 long). I tried to limit all those potential problems by putting additionnal assertions on the sequence, the strand and the position we ask (depending of course of the indexing system). All of these security assumptions are implemented. I just havn't really had the time to test it yet, so be carefull with that.

#  <a name="whichmodel"></a> Which model should I use ?

Available models are:

* C2
* P4-C2
* P5-C3
* P6-C4
* SP2-C2
* XL-C2
* XL-XL

**No model** is available for **Sequel II** with **chemistry V1.0**

**SP2-C2 is the latest** of all models [Jan 2020], recommended for **Sequel I and Sequel II v2 data**, including for people using P6-C4 chemistry

**Most RSII user would probably want to use P5-C3 model**

> See [here](https://github.com/PacificBiosciences/kineticsTools/pull/71) for more info.


# Conventions used

> **Important notes:**
> - 0-based or 1-based indexing depends on the user's choice (PacBio uses 1-based so, our default is 0-based unless specified explicitely)
> - Position N strand 0 refers to the Nth position in the fasta
> - **Strand 0 is the fasta reference sequence**. Strand 1 is the complementary
> - Position is **alswyas** expressed **from the strand0/fasta reference point of view**. We don't care about 5'-3'.

E.g, if we use 0-based indexing (not the default one) then,

Given this fasta sequence :

```console
"ATGCTT"
```
The DNA molecule could be pictured this way :

<pre>
position                               012345
strand0                             5' ATGCTT 3'
                                       ||||||
strand1                             3' TACGAA 5'
position                               012345
</pre>

Which implies, for instance that:

- Prediction of position=2 and strand = 1 will be the prediction at position 2 on the complementary strand: "C"
- Prediction of position=2 and strand = 0 will be the prediction at position 2 on read strand : "G"

**Here is a table for the output that will be produced from this exemple just to be sure:**

For strand 0

position | strand | nucleotide |
:-:|:-:|:-:|
0|0|A
1|0|T
2|0|G
3|0|C
4|0|T
5|0|T

And for strand 1

position | strand | nucleotide |
:-:|:-:|:-:|
0|1|T
1|1|A
2|1|C
3|1|G
4|1|A
5|1|A


# Credits

The Original code is copyrighted from Pacific Bioscience California (See "LICENSE") as it is published on Github:  [kineticsTools](https://github.com/PacificBiosciences/kineticsTools) at the date of 30 January 2020.

Major modifications, retro-ingineering, isolation from the PacBio's API, repackaging, python3 compatibility and new documentation are from [Guillaume DELEVOYE](https://github.com/GDelevoye).

See "License" for more infos
