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

> Works well on Ubuntu 18.04 LTS x86. Should'nt work on MacOS or Windows (A dynamical library will be missing)

# Installation

ipdtools can be installed with **pip**, which is automatically distributed with any version of python

```bash
git clone https://github.com/GDelevoye/ipdtools.git
pip install ./ipdtools
```

Please note that, due to the fact we're compiling a third-party dependency, the "editable" pip install will not be available. Do not use pip install -e ./ipdtools

# Usage

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
                [--verbosity {DEBUG,INFO,WARNING,ERROR, CRITICAL}]

optional arguments:
  -h, --help            show this help message and exit
  --model {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL}, -m {SP2-C2,C2,P4-C2,P5-C3,P6-C4,XL-C2,XL-XL}
                        Choose the model for IPD prediction. See the README of
                        package for more info. DEFAULT: SP2-C2
  --fastafile FASTAFILE, -f FASTAFILE
                        Path to a fasta file.
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        Output CSV file of predicted IPDs.
  --verbosity {DEBUG,INFO,WARNING,ERROR, CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR, CRITICAL}
```

> RQ: If verbosity is DEBUG or INFO, the program will show a tqdm progress bar

```console

```(py2) user@computer:~/ipdtools/ipdtools/resources$ ipdtools -f test.fasta -o test.csv -v DEBUG
100%|██████████| 93/93 [00:00<00:00, 2954.90it/s]
100%|██████████| 92/92 [00:00<00:00, 2890.54it/s]
100%|██████████| 92/92 [00:00<00:00, 2919.34it/s]
100%|██████████| 93/93 [00:00<00:00, 2937.36it/s]
100%|██████████| 93/93 [00:00<00:00, 2935.88it/s]
```

### Output .csv header

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

> The csv is sorted (in order) by:
> * Fasta_ID
> * Position
> * Strand

## Python API

ipdtools can be used in python to predict from a fasta file


```python
import ipdtools

fasta = ipdtools.ipdModel.loadReferenceContigs('test.fasta')
model = ipdtools.ipdModel.IpdModel(fasta,modelname="SP2-C2")

position = 15
strand = 0

prediction = model.predictIpdFunc(identifier="seq0")(position,strand)
print(prediction)
```

it can also be used to predict directly from a python string with **Str2IPD**:

```python
>>> import ipdtools
>>> predictor = ipdtools.ipdModel.Str2IPD("ATGCTAGCTTTTTGNCTGATTAGCTGA",model="SP2-C2") # Default model is SP2-C2
>>> predictor.predict(position=0) # The default is strand0
1.0641327
>>> predictor.predict(position=0,strand=1)
0.42285475
```

Because most of the times you want a prediction for only one model, but many sequences, you can use **batchStr2IPD** which works on the same principle. Here is it's prototype:

```python
batchStr2IPD(
              sequences=["ATCGATGCTGGCGCT",
                         "ATCGTAGCTGATCGGAGTG",
                         "CGCGCGCGCGTTA"],
              names=["seq1","seq2","seq3"],
              model="SP2-C2")
```

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
> - We are using **0-indexing**
> - Position N strand 0 refers to the Nth position in the fasta
> - **Strand 0 is the fasta file**. Strand 1 is the complementary
> - Position is **alswyas** expressed **from the strand0 point of view**. We don't care about 5'-3'.

E.g

Given this fasta sequence :

```console
"ATGCTT"
```
We then have in real :

<pre>
position                               012345
strand0                             5' ATGCTT 3'
                                       ||||||
strand1                             3' TACGAA 5'
position                               012345
</pre>

e.g

- Prediction of position=2 and strand = 1 will be the prediction at position 2 on the complementary strand: "C"
- Prediction of position=2 and strand = 0 will be the prediction at position 2 on read strand : "G"

**Here is a table for this exemple just to be sure:**

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
