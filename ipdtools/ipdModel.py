#!/usr/bin/env python
################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
################################################################################

#### Major modifications and monkeypatching operated by DELEVOYE Guillaume


import logging
import os
import h5py
import ctypes as C
from pkg_resources import Requirement, resource_filename

from multiprocessing.sharedctypes import RawArray
import warnings
import numpy as np

import pandas as pd
from joblib import Parallel, delayed
import multiprocessing as mp
from tqdm import tqdm

import gc # Garbage Collector
import sys

byte = np.dtype('byte')
float32 = np.dtype('float32')
uint8 = np.dtype('uint8')

# Map for ascii encoded bases to integers 0-3 -- will be used to define a 24-bit lookup code
# for fetching predicted IPDs from the kinetic LUT.

# We start everything at 0, so anything will map to 'A' unless it appears in this table
lutCodeMap = np.zeros(256, dtype=uint8)
maps = {'a': 0, 'A': 0, 'c': 1, 'C': 1, 'g': 2, 'G': 2, 't': 3, 'T': 3}
for k in maps:
    lutCodeMap[ord(k)] = maps[k]
lutReverseMap = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

seqCodeMap = np.ones(256, dtype=uint8) * 4
for k in maps:
    seqCodeMap[ord(k)] = maps[k]
seqMap = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}
seqMapNp = np.array(['A', 'C', 'G', 'T', 'N'])

seqMapComplement = {0: 'T', 1: 'G', 2: 'C', 3: 'A', 4: 'N'}
seqMapComplementNp = np.array(['T', 'G', 'C', 'A', 'N'])

# Base letters for modification calling
# 'H' : m6A, 'I' : m5C, 'J' : m4C, 'K' : m5C/TET
baseToCode = {'N': 0, 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'H': 4, 'I': 5, 'J': 6, 'K': 7}
baseToCanonicalCode = {'N': 0, 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'H': 0, 'I': 1, 'J': 1, 'K': 1}

codeToBase = dict([(y, x) for (x, y) in list(baseToCode.items())])

def result_writer(queue, fileOut, length, show_progress_bar):
    """ Will be launched as an independant process, writing every that will work until it reads 'poison' on the queue. After it's done, he will close the files"""

    os.system('mkdir -p '+str(os.path.dirname(fileOut)))
    list_df = []

    logging.debug('[DEBUG] (result_writer) result_writer is ready to write (destination = {})'.format(fileOut))
    logging.info("[INFO] Number of nucleotides to handle : {}".format(str(length)))


    if show_progress_bar:
        progress_bar = tqdm(range(length),total=length,desc="Nucleotides")


    while True:  # Continuously listens to the queue
        if not queue.empty():  # most of the time it will be empty
            to_write = queue.get() # This is a pandas DataFrame otherwise the word "poison" when we reached the end of analysis

            if not isinstance(to_write,pd.DataFrame):  # If the queue is poisoned
                final_result = pd.concat(list_df, ignore_index = True)

                with open(os.path.realpath(fileOut),"w") as myfile:
                    final_result.to_csv(myfile, sep = ';', index = False)
                    #logging.debug('[DEBUG] (result_writer) Written in {}'.format(fileOut))

                break

            else:
                list_df.append(to_write.copy())
                if show_progress_bar:
                        progress_bar.update(len(to_write)//2) # Don't forget //2 because we have two strands

    sys.stdout.flush()
    logging.info('[INFO] Result saved at {}'.format(os.path.realpath(fileOut)))
    sys.stdout.flush()

    return

def dict_to_fasta(myDict, fastafile):
    """ from a python dict {Identifier:sequence}, writes a fastafile"""
    def divide_seq(seq, length = 60, sep = '\n'):
        """From a full-length sequence, divides it in chunks of 60 letters (most frequent.fasta format)"""
        return( str(sep).join([seq[x:x+int(length)] for x in range(0, len(seq), length)]))
        # Tricky one liner --> Splits a sequence into chunks of size "length"
    with open(os.path.realpath(fastafile), "w") as filout:
        for key in myDict:
            filout.write('>'+str(key)+'\n'+str(divide_seq(myDict[key]))+'\n')


def parallelized_function(contig_name,fastafile,modelname,queue,indexing):

    logging.debug('[DEBUG] Caring about contig : {}'.format(contig_name))
    fasta = load_fasta(os.path.realpath(fastafile)) # Same
    seq = fasta[contig_name]

    # Creating a temporary fasta file (A bit dirty, but way quicker than having to rework the code completely)
    tmp_fasta = "./tmp"+contig_name+".fasta"
    dict_to_fasta({contig_name:seq}, tmp_fasta)
    fastafile = tmp_fasta

    
    logging.debug("[DEBUG] contig {}, Creating a model from FastaRecords".format(contig_name))
    fastaRecords = loadReferenceContigs(os.path.realpath(fastafile)) ### Not optimal to do it each time here
    # But the object "model" cannot be pickled so it has to be done this way
    cpl = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    logging.debug("[DEBUG] contig {}, Creating a model from FastaRecords".format(contig_name))
    model = IpdModel(fastaRecords,modelFile=transform_model_name(modelname))

    logging.debug("[DEBUG] contig {}, Fetching the sequence".format(contig_name))

    logging.debug("[DEBUG] Caring about contig {} of size {}".format(contig_name,len(seq)))

    predictfunc = model.predictIpdFuncModel(refId=contig_name)

    list_return = []

    for i in range(len(seq)):  # Progress bar will be shown only if specified
        prediction_strand0 = float(predictfunc(i, 0))  # Prise en compte de l'indexing 1-based par PacBio
        prediction_strand1 = float(predictfunc(i, 1))

        if indexing == 1:
            list_return.append({"Fasta_ID": contig_name, "Position": i+1, "Strand": 0, "Nucleotide": str(seq[i]),
                                "Prediction": float(prediction_strand0)})
            list_return.append({"Fasta_ID": contig_name, "Position": i+1, "Strand": 1, "Nucleotide": str(cpl[seq[i]]),
                                "Prediction": float(prediction_strand1)})
        elif indexing == 0:

            list_return.append({"Fasta_ID": contig_name, "Position": i, "Strand": 0, "Nucleotide": str(seq[i]),
                                "Prediction": float(prediction_strand0)})
            list_return.append({"Fasta_ID": contig_name, "Position": i, "Strand": 1, "Nucleotide": str(cpl[seq[i]]),
                                "Prediction": float(prediction_strand1)})

    queue.put(pd.DataFrame(list_return).sort_values(["Fasta_ID","Position","Strand"]).copy())

    os.remove(os.path.realpath(tmp_fasta))
    del list_return
    gc.collect()

def compute_fasta_to_csv(modelname, fastafile, csv_out, show_progress_bar=False, nproc=1,indexing=1):
    fasta = load_fasta(os.path.realpath(fastafile))

    length = sum([len(fasta[x]) for x in fasta])

    contig_names = [contig_name for contig_name in fasta]

    logging.info("[INFO] Using {} CPUs".format(nproc))

    logging.debug('[DEBUG] (compute_fasta_to_csv) Creating the manager queue')
    manager = mp.Manager()
    queue = manager.Queue()

    logging.debug('[DEBUG] (compute_fasta_to_csv) Launching the writer_process')
    writer_process = mp.Process(target=result_writer, args=(queue,os.path.realpath(csv_out),length,show_progress_bar))
    writer_process.start()  # This will write every worker's result that has been pushed in the queue

    logging.debug('[DEBUG] (compute_fasta_to_csv) Waiting for the workers to be done')


    Parallel(n_jobs=nproc)(
        delayed(parallelized_function)(contig_name, fastafile, modelname, queue, indexing) for contig_name in contig_names)

    ################################################
    # END OF THE ANALYSIS / WAITING FOR THE WRITER #
    ################################################

    queue.put("Poison")  # Poisoning the end of the queue
    writer_process.join()  # Waiting writer_process to reach the poison
    #logging.info('[INFO] (compute_fasta_to_csv) Analysis DONE and saved')  # You're welcome

    logging.info('[INFO] DONE')

class Str2IPD():
    def __init__(self, sequence, name="NO_ID", model="SP2-C2",indexing=1):
        assert not any([character not in ["A","T","C","G","N","a","t","c","g","n"] for character in sequence])
        assert indexing in [0,1]

        self.indexing = indexing
        assert sequence != None and sequence != "" and len(sequence) > 0
        self.character_seq = sequence

        self.sequence = [Contig(name, sequence)]
        for x in self.sequence:
            x.cmph5ID = x.id

        self.model = IpdModel(self.sequence, modelFile=transform_model_name(model))
        self.predictfunc = self.model.predictIpdFuncModel(refId=name)

    def predict(self, position, strand=0):
        assert isinstance(position,int)
        assert strand in [0,1]

        if self.indexing == 1:
            assert position >= 1
            assert position+1 <= len(self.character_seq)
            return self.predictfunc(position+1, strand)
        elif self.indexing == 0:
            assert position >= 0
            assert position <= len(self.character_seq)
            return self.predictfunc(position, strand)



class batchStr2IPD():
    def __init__(self, sequences=[], names=[], model="SP2-C2"):
        self.sequences = [Contig(name, sequence) for (name, sequence) in zip(sequences, names)]
        for x in self.sequences:
            x.cmph5ID = x.id

        self.model = IpdModel(self.sequences, modelFile=transform_model_name(model))
        self.predictfunc = {x: self.model.predictIpdFuncModel(refId=name) for name in names}

    def predict(self, sequence, position, strand=0):
        return self.predictfunc[sequence](position, strand)


def transform_model_name(modelname):
    resources_dir = _getAbsPath("/resources/")
    modifiedmodelname = modelname + ".h5"
    modelname = os.path.join(resources_dir, modifiedmodelname)
    return modelname


def load_fasta(fastafile):
    """Returns a python dict { id : sequence } for the given .fasta file"""
    with open(os.path.realpath(fastafile), 'r') as filin:
        fasta = filin.read()
        fasta = fasta.split('>')[1:]
        outputdict = {x.split('\n')[0].strip(): "".join(x.split('\n')[1:]) for x in fasta}
    return outputdict


class Contig():
    def __init__(self, id, seq):
        self.id = id
        self.ID = id
        self.sequence = seq


def loadReferenceContigs(referencePath):
    """
    Load the reference contigs, and tag each one with the ref.cmpH5ID it
    was assigned in the alignment file(s).  Return a list of contigs,
    which are used to set up IpdModel.
    """

    fasta_parsed = load_fasta(os.path.realpath(referencePath))
    contigs = [Contig(x, fasta_parsed[x]) for x in fasta_parsed.keys()]

    # Read contigs from FASTA file (or XML dataset)
    # refReader = ReferenceSet(referencePath)
    # contigs = []
    # contigs.extend([x for x in refReader])
    # contigDict = dict([(x.id,  x) for x in contigs])

    # initially each contig has an id of None -- this will be overwritten with the id from the cmp.h5, if there are any
    # reads mapped to it.
    for x in contigs:
        x.cmph5ID = x.id

    # # Mark each contig with it's ID from the cmp.h5 - match them up using MD5s
    # for x in alignmentSet.referenceInfoTable:
    #     if x.FullName in contigDict:
    #         contigDict[x.FullName].cmph5ID = x.ID

    return contigs


class SharedArray:
    """
    Very simple wrapper for a chunk of shared memory that can be accessed across processes
    """

    def __init__(self, dtype, shape):
        self._rawArray = RawArray(dtype, shape)

    def getNumpyWrapper(self):
        """
        Construct a numpy array that wraps the raw shared memory array
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.ctypeslib.as_array(self._rawArray)


def _getAbsPath(fname):
    return resource_filename(Requirement.parse('ipdtools'), 'ipdtools/%s' % fname)


class GbmContextModel(object):
    """
    Class for computing ipd predictions on contexts. Evaluate the GBM tree model for a list of contexts
    Contexts may contain arbitrary combinations of modified bases
    """

    def __init__(self, modelH5Group, modelIterations=-1):

        # This will hold the ctypes function pointer
        # It will be lazily initialized
        self.nativeInnerPredict = None
        self.nativeInnerPredictCtx = None

        def ds(name):
            return modelH5Group[name][:]

        # print("ds(Varnames)",ds("VarNames"))
        self.varNames = [x.decode('ASCII') for x in ds("VarNames")]
        # print("self.varNames2",self.varNames)

        self.modFeatureIdx = dict(
            (int(self.varNames[x][1:]), x) for x in range(len(self.varNames)) if self.varNames[x][0] == 'M')
        self.canonicalFeatureIdx = dict(
            (int(self.varNames[x][1:]), x) for x in range(len(self.varNames)) if self.varNames[x][0] == 'R')

        # print("self.modFeatureIdx: {}".format(self.modFeatureIdx))
        # print("self.canonicalFeatureIdx",self.canonicalFeatureIdx)

        self.pre = 10
        self.post = 4
        self.ctxSize = self.pre + self.post + 1

        self.splitVar = ds("Variables")
        self.leftNodes = ds("LeftNodes")
        self.rightNodes = ds("RightNodes")
        self.missingNodes = ds("MissingNodes")

        self.splitVar16 = self.splitVar.astype(np.int16)

        self.splitCodes = ds("SplitCodes").astype(np.float32)

        self.cSplits = ds("CSplits")
        self.maxCSplits = self.cSplits.shape[1]

        self.initialValue = ds("InitialValue").astype(np.float32)[0]

        exp = 2 ** np.arange(self.cSplits.shape[1] - 1, -1, -1)
        self.bSplits = ((self.cSplits > 0) * exp).sum(1)

        # total number of trees in model
        self.nTrees = self.splitVar.shape[0]
        self.treeSize = self.splitVar.shape[1]

        offsets = np.floor(np.arange(0, self.leftNodes.size) / self.treeSize) * self.treeSize
        offsets = offsets.astype(np.int32)

        self.leftNodesOffset = self.leftNodes.flatten().astype(np.int32) + offsets
        self.rightNodesOffset = self.rightNodes.flatten().astype(np.int32) + offsets
        self.missingNodesOffset = self.missingNodes.flatten().astype(np.int32) + offsets

        self.splitCodesCtx = self.splitCodes.copy().flatten()

        splitCodesCtxView = self.splitCodesCtx.view()
        splitCodesCtxView.dtype = np.uint32

        # Pack the cSplits as a bit array directly into the splitCode array
        # using an uin32 view of the splitCode array
        flatSplitVar = self.splitVar.flatten()

        powOfTwo = 2 ** np.arange(self.maxCSplits)

        for i in range(self.splitCodesCtx.shape[0]):

            if flatSplitVar[i] != -1:
                # This is a pointer to a cSplit row -- pack the csplit into a unit32, then overwirte
                # this slot of the ctxSplitCodes
                cs = self.cSplits[int(self.splitCodesCtx[i]), :]
                v = (powOfTwo * (cs > 0)).sum()

                splitCodesCtxView[i] = v

        # If the user has requested fewer iterations, update nTrees
        if modelIterations > 0:
            self.nTrees = modelIterations

    def _initNativeTreePredict(self):
        """
        Initialization routine the C tree-predict method
        Needs to be invoked lazily because the native function pointer cannot be pickled
        """

        DLL_PATH = os.path.dirname(os.path.dirname(__file__) + "/ipdtools/ipdtools")
        # print("******************"+DLL_PATH+"******************")

        self._lib = np.ctypeslib.load_library("tree_predict", DLL_PATH)

        lpb = self._lib

        lpb.init_native.argtypes = [C.c_int]

        fp = C.POINTER(C.c_float)
        fpp = C.POINTER(fp)
        ip = C.POINTER(C.c_int)
        sp = C.POINTER(C.c_int16)
        ui64p = C.POINTER(C.c_uint64)

        args = [fp, fpp, C.c_int, ip, ip, ip, fp, ip, ip, ip, C.c_float, C.c_int, C.c_int, C.c_int]
        lpb.innerPredict.argtypes = args
        self.nativeInnerPredict = lpb.innerPredict

        # Fast version

        # void innerPredictCtx(
        #    int ctxSize, float radPredF[], uint64_t contextPack[], int cRows,
        #    int16 left[], int16 right[], int16 missing[], float splitCode[], int16 splitVar[],
        #    int varTypes[], float initialValue, int treeSize, int numTrees, int maxCSplitSize)

        args = [C.c_int, fp, ui64p, C.c_int, ip, ip, ip, fp, sp, ip, C.c_float, C.c_int, C.c_int, C.c_int]
        lpb.innerPredictCtx.argtypes = args
        self.nativeInnerPredictCtx = lpb.innerPredictCtx

    def getPredictionsSlow(self, ctxStrings, nTrees=None):
        """Compute IPD predictions for arbitrary methylation-containing contexts."""
        # C prototype that we call:
        # void innerPredict(
        #   float[] radPredF,
        #   IntPtr[] dataMatrix,
        #   int cRows, int[] left, int[] right, int[] missing,
        #   float[] splitCode, int[] splitVar, int[] cSplits,
        #   int[] varTypes, float initialValue,
        #   int treeSize, int numTrees, int maxCSplitSize);

        # Make sure native library is initialized
        if self.nativeInnerPredict is None:
            self._initNativeTreePredict()

        def fp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_float))

        def ip(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int))

        if nTrees is None:
            nTrees = self.nTrees

        n = len(ctxStrings)

        mCols = [np.zeros(n, dtype=np.float32) for x in range(self.ctxSize)]
        rCols = [np.zeros(n, dtype=np.float32) for x in range(self.ctxSize)]

        for stringIdx in range(len(ctxStrings)):
            s = ctxStrings[stringIdx]

            for i in range(len(s)):
                mCols[i][stringIdx] = baseToCode[s[i]]
                rCols[i][stringIdx] = baseToCanonicalCode[s[i]]

        dataPtrs = (C.POINTER(C.c_float) * (2 * self.ctxSize))()

        varTypes = np.zeros(2 * self.ctxSize, dtype=np.int32)

        for i in range(self.ctxSize):
            dataPtrs[self.modFeatureIdx[i]] = mCols[i].ctypes.data_as(C.POINTER(C.c_float))
            dataPtrs[self.canonicalFeatureIdx[i]] = rCols[i].ctypes.data_as(C.POINTER(C.c_float))

            varTypes[self.modFeatureIdx[i]] = 8
            varTypes[self.canonicalFeatureIdx[i]] = 4

        self.predictions = np.zeros(len(ctxStrings), dtype=np.float32)

        self.nativeInnerPredict(
            fp(self.predictions), dataPtrs,
            n, ip(self.leftNodes), ip(self.rightNodes), ip(self.missingNodes),
            fp(self.splitCodes), ip(self.splitVar), ip(self.cSplits),
            ip(varTypes), self.initialValue, self.treeSize, nTrees, self.maxCSplits)

        return np.exp(self.predictions)

    def getPredictions(self, ctxStrings, nTrees=None):
        """Compute IPD predictions for arbitrary methylation-containing contexts."""
        # C prototype that we call:
        # void innerPredictCtx(
        #   int ctxSize, float[] radPredF,
        #   int[] contextPack,
        #   int cRows, int[] left, int[] right, int[] missing,
        #   float[] splitCode, int[] splitVar,
        #   int[] varTypes, float initialValue,
        #   int treeSize, int numTrees, int maxCSplitSize);

        # Make sure native library is initialized
        if self.nativeInnerPredictCtx is None:
            self._initNativeTreePredict()

        def fp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_float))

        def ip(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int))

        def ulp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_uint64))

        def sp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int16))

        n = len(ctxStrings)

        if nTrees is None:
            nTrees = self.nTrees

        packCol = np.zeros(n, dtype=np.uint64)

        for stringIdx in range(len(ctxStrings)):
            s = ctxStrings[stringIdx].decode('UTF-32')
            code = 0

            for i in range(len(s)):
                # print(s.decode("utf-8"),i,s.decode('utf-8')[i],baseToCode) # For debug purposes
                # print("baseToCode",baseToCode)
                # print("s",s)
                # print("s",s,"i",i,"s[i]",s[i],"baseToCode",baseToCode,"code",code)
                # print(s[i])
                # print("".join([chr(x) for x in s]),chr(s[i]))
                # print(s,i,s[i])
                modBits = baseToCode[s[i]]

                # print(i,self.modFeatureIdx)
                slotForPosition = self.modFeatureIdx[i]
                # print("modBits",modBits,"slotForPosition",slotForPosition)

                code = code | (modBits << (4 * slotForPosition))

            packCol[stringIdx] = code

        # print packed base codes
        # for v in packCol.flatten():
        #    print v
        #    for i in np.arange(12):
        #        print "%d: %o" % (i,  (v.item() >> (5*i)) & 0x1f)

        varTypes = np.zeros(2 * self.ctxSize, dtype=np.int32)

        for i in range(self.ctxSize):
            varTypes[self.modFeatureIdx[i]] = 8
            varTypes[self.canonicalFeatureIdx[i]] = 4

        self.predictions = np.zeros(len(ctxStrings), dtype=np.float32)

        self.nativeInnerPredictCtx(
            self.ctxSize, fp(self.predictions), ulp(packCol),
            n, ip(self.leftNodesOffset), ip(self.rightNodesOffset), ip(self.missingNodesOffset),
            fp(self.splitCodesCtx), sp(self.splitVar16),
            ip(varTypes), self.initialValue, self.treeSize, nTrees, self.maxCSplits)

        return np.exp(self.predictions)


class IpdModel:
    """
    Predicts the IPD of an any context, possibly containing multiple modifications.
    We use a 4^12 entry LUT to get the predictions for contexts without modifications,
    then we use the GbmModel to get predictions in the presence of arbitrary mods.
        Note on the coding scheme.  For each contig we store a byte-array that has size = contig.length + 2*self.pad
        The upper 4 bits contain a lookup into seqReverseMap, which can contains N's. This is used for giving
        template snippets that may contains N's if the reference sequence does, or if the snippet
        The lowe 4 bits contain a lookup into lutReverseMap, which
    """

    def __init__(self, fastaRecords, modelFile, modelIterations=-1):
        """
        Load the reference sequences and the ipd lut into shared arrays that can be
        used as numpy arrays in worker processes.
        fastaRecords is a list of FastaRecords, in the cmp.h5 file order
        """

        self.pre = 10
        self.post = 4

        self.pad = 30
        self.base4 = 4 ** np.array(list(range(self.pre + self.post + 1)))

        self.refDict = {}
        self.refLengthDict = {}

        for contig in fastaRecords:
            rawSeq = contig.sequence[:]
            refSeq = np.fromstring(rawSeq, dtype=byte)

            # Store the reference length
            self.refLengthDict[contig] = len(rawSeq)

            # Make a shared array
            sa = SharedArray(dtype='B', shape=len(rawSeq) + self.pad * 2)
            saWrap = sa.getNumpyWrapper()

            # Lut Codes convert Ns to As so that we don't put Ns into the Gbm Model
            # Seq Codes leaves Ns as Ns for getting reference snippets out
            innerLutCodes = lutCodeMap[refSeq]
            innerSeqCodes = seqCodeMap[refSeq]
            innerCodes = np.bitwise_or(innerLutCodes, np.left_shift(innerSeqCodes, 4))

            saWrap[self.pad:(len(rawSeq) + self.pad)] = innerCodes

            # Padding codes -- the lut array is padded with 0s the sequence array is padded with N's (4)
            outerCodes = np.left_shift(np.ones(self.pad, dtype=uint8) * 4, 4)
            saWrap[0:self.pad] = outerCodes
            saWrap[(len(rawSeq) + self.pad):(len(rawSeq) + 2 * self.pad)] = outerCodes

            self.refDict[contig.cmph5ID] = sa

        # No correction factor for IPDs everything is normalized to 1
        self.meanIpd = 1

        # Find and open the ipd model file
        self.lutPath = modelFile
        if os.path.exists(self.lutPath):
            h5File = h5py.File(self.lutPath, mode='r')

            gbmModelGroup = h5File["/AllMods_GbmModel"]
            self.gbmModel = GbmContextModel(gbmModelGroup, modelIterations)

            # We always use the model -- no more LUTS
            self.predictIpdFunc = self.predictIpdFuncModel
            self.predictManyIpdFunc = self.predictManyIpdFuncModel
        else:
            logging.info("Couldn't find model file: %s" % self.lutPath)

    def _loadIpdTable(self, nullModelGroup):
        """
        Read the null kinetic model into a shared numpy array dataset
        """
        nullModelDataset = nullModelGroup["KineticValues"]

        # assert that the dataset is a uint8
        assert (nullModelDataset.dtype == uint8)

        # Construct a 'shared array' (a numpy wrapper around some shared memory
        # Read the LUT into this table
        self.sharedArray = SharedArray('B', nullModelDataset.shape[0])
        lutArray = self.sharedArray.getNumpyWrapper()
        nullModelDataset.read_direct(lutArray)

        # Load the second-level LUT
        self.floatLut = nullModelGroup["Lut"][:]

    def refLength(self, refId):
        return self.refLengthDict[refId]

    def cognateBaseFunc(self, refId):
        """
        Return a function that returns a snippet of the reference sequence around a given position
        """

        # FIXME -- what is the correct strand to return?!
        # FIXME -- what to do about padding when the snippet runs off the end of the reference
        # how do we account for / indicate what is happening
        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, tplStrand):

            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand == 0:
                slc = refArray[tplPos]
                slc = np.right_shift(slc, 4)
                return seqMap[slc]

            # Reverse strand
            else:
                slc = refArray[tplPos]
                slc = np.right_shift(slc, 4)
                return seqMapComplement[slc]

        return f

    def snippetFunc(self, refId, pre, post):
        """
        Return a function that returns a snippet of the reference sequence around a given position
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, tplStrand):
            """Closure for returning a reference snippet. The reference is padded with N's for bases falling outside the extents of the reference"""
            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand == 0:
                slc = refArray[(tplPos - pre):(tplPos + 1 + post)]
                slc = np.right_shift(slc, 4)
                return seqMapNp[slc].tostring()

            # Reverse strand
            else:
                slc = refArray[(tplPos + pre):(tplPos - post - 1):-1]
                slc = np.right_shift(slc, 4)
                return seqMapComplementNp[slc].tostring()

        return f

    def getReferenceWindow(self, refId, tplStrand, start, end):
        """
        Return  a snippet of the reference sequence
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        # adjust position for reference padding
        start += self.pad
        end += self.pad

        # Forward strand
        if tplStrand == 0:
            slc = refArray[start:end]
            slc = np.right_shift(slc, 4)
            return "".join(seqMap[x] for x in slc)

        # Reverse strand
        else:
            slc = refArray[end:start:-1]
            slc = np.right_shift(slc, 4)
            return "".join(seqMapComplement[x] for x in slc)

    def predictIpdFuncLut(self, refId):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        # Materialized the numpy wrapper around the shared data
        refArray = self.refDict[refId].getNumpyWrapper()
        lutArray = self.sharedArray.getNumpyWrapper()
        floatLut = self.floatLut

        def f(tplPos, tplStrand):

            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand == 0:
                slc = np.bitwise_and(refArray[(tplPos + self.pre):(tplPos - self.post - 1):-1], 0xf)

            # Reverse strand
            else:
                slc = 3 - np.bitwise_and(refArray[(tplPos - self.pre):(tplPos + 1 + self.post)], 0xf)

            code = (self.base4 * slc).sum()
            return floatLut[max(1, lutArray[code])]

        return f

    def predictIpdFuncModel(self, refId):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        # Materialized the numpy wrapper around the shared data
        snipFunction = self.snippetFunc(refId, self.post, self.pre)

        def f(tplPos, tplStrand):
            # Get context string
            context = snipFunction(tplPos, tplStrand)

            # Get prediction
            return self.gbmModel.getPredictions([context])[0]

        return f

    def predictManyIpdFuncModel(self, refId):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        # Materialized the numpy wrapper around the shared data
        snipFunction = self.snippetFunc(refId, self.post, self.pre)

        def fMany(sites):
            contexts = [snipFunction(x[0], x[1]) for x in sites]
            return self.gbmModel.getPredictions(contexts)

        return fMany

    def modPredictIpdFunc(self, refId, mod):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, relativeModPos, readStrand):

            # skip over the padding
            tplPos += self.pad

            # Read sequence matches forward strand
            if readStrand == 0:
                slc = 3 - np.bitwise_and(refArray[(tplPos - self.pre):(tplPos + 1 + self.post)], 0xf)

            # Reverse strand
            else:
                slc = np.bitwise_and(refArray[(tplPos + self.pre):(tplPos - self.post - 1):-1], 0xf)

            # Modify the indicated position
            slc[relativeModPos + self.pre] = baseToCode[mod]

            slcString = "".join([codeToBase[x] for x in slc])

            # Get the prediction for this context
            # return self.gbmModel.getPredictions([slcString])[0]
            return 0.0

        return f
