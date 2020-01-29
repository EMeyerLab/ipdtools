# Warning

Due to encoding, only work in python2

# ipdtools
Retro-ingeniering of PacBio's modelling of IPD (RSII, Sequel I, Sequel II v2)

```console
(py2) delevoye@gmcpc04:~/Github/ipdtools/ipdtools/resources$ python
Python 2.7.17 |Anaconda, Inc.| (default, Oct 21 2019, 19:04:46) 
[GCC 7.3.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
```
```python
>>> import ipdtools
>>> fastaRecords= ipdtools.ipdModel.loadReferenceContigs('test.fasta')
>>> model = ipdtools.ipdModel.IpdModel(fastaRecords,"./SP2-C2.h5")
>>> model.predictIpdFunc("seq5")(15,0)
1.4750001
>>> model.predictIpdFunc("seq5")(0,0)
1.3652387
>>>
```

