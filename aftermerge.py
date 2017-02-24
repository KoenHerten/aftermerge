#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 14:03:19 2017

@author: Koen Herten 

This is aftermerge v1. A script for merging paired reads after mapping, based on the locations on the reference genome.

Copyright 2017, Koen Herten, All rights reserved

This file is part of aftermerge.

aftermerge is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

aftermerge is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with aftermerge.  If not, see <http://www.gnu.org/licenses/>.
"""

from argparse import ArgumentParser
from argparse import FileType
import sys
import samRead
   

def mergeSequences(seq1, seq2, qual1, qual2, mcigar1, mcigar2, start1, start2):
    '''
    function returns (sequence, quality) the merged sequence and quality
    when bases in the merge match, base quality is the max
    when they do not match, base is base with highest quality, quality is the average
    '''
    seq1 = str(seq1)
    seq2 = str(seq2)
    qual1 = str(qual1)
    qual2 = str(qual2)
    mcigar1 = str(mcigar1)
    mcigar2 = str(mcigar2)
    start1 = int(start1)
    start2 = int(start2)
    newseq = ""
    newqual = ""
    newmcigar = ""
    if (seq1 == ""):
        #no seq1 any more
        return (seq2, qual2, mcigar2)
    if (seq2 == ""):
        #no seq2 any more
        return (seq1, qual1, mcigar1)
    if (start2 < start1):
        #seq2 is before seq1
        (newseq, newqual, newmcigar) = mergeSequences(seq1, seq2[1:], qual1, qual2[1:], mcigar1, mcigar2[1:], start1, start2+1)
        newseq = seq2[0] + newseq
        newqual = qual2[0] + newqual
        newmcigar = mcigar2[0] + newmcigar
    elif (start2 > start1):
        #seq1 is before seq2
        (newseq, newqual, newmcigar) = mergeSequences(seq1[1:], seq2, qual1[1:], qual2, mcigar1[1:], mcigar2, start1+1, start2)
        newseq = seq1[0] + newseq
        newqual = qual1[0] + newqual
        newmcigar = mcigar1[0] + newmcigar
    else:
        #same pos
        if (mcigar1[0] == mcigar2[0]):
            #same cigar, no problem
            (newseq, newqual, newmcigar) = mergeSequences(seq1[1:], seq2[1:], qual1[1:], qual2[1:], mcigar1[1:], mcigar2[1:], start1+1, start2+1)
            if (seq1[0] == seq2[0]):
                #same base, so take highest quality
                newseq = seq1[0] + newseq
                qual = max(qual1[0], qual2[0])
                newqual = qual + newqual
                newmcigar = mcigar1[0] + newmcigar
            else:
                #not same base, so take base with highest quality, give lowest quality
                base = seq1[0]
                if (qual1[0] < qual2[0]):
                    base = seq2[0]
                qual = min(qual1[0], qual2[0])
                newseq = base + newseq
                newqual = qual + newqual
                newmcigar = mcigar1[0] + newmcigar
        else:
            #not same cigar
            c1 = str(mcigar1[0])
            c2 = str(mcigar2[0])
            if (c1.isupper() and c2.isupper()):
                #both match or mismatch
                (newseq, newqual, newmcigar) = mergeSequences(seq1[1:], seq2[1:], qual1[1:], qual2[1:], mcigar1[1:], mcigar2[1:], start1+1, start2+1)
                base = seq1[0]
                m = mcigar1[0]
                if (qual1[0] < qual2[0]):
                    base = seq2[0]
                    m = mcigar2[0]
                qual = chr(int((ord(qual1[0]) + ord(qual2[0])) / 2))
                newseq = seq1[0] + newseq
                newqual = qual + newqual
                newmcigar = m + newmcigar
            else:
                #one is insert, deletion or clipped
                if (c1 == "s"):
                    #seq1 is softclipped
                    (newseq, newqual, newmcigar) = mergeSequences("", seq2, "", qual2, "", mcigar2, start1+1, start2+1)
                elif (c2 == "s"):
                    #seq2 is softclipped
                    (newseq, newqual, newmcigar) = mergeSequences(seq1, "", qual1, "", mcigar1, "", start1+1, start2+1)
                else:
                    #One is insertion or deletion
                    print("unable to resolve: {}\t{}".format(c1, c2))
    return (newseq, newqual, newmcigar)

    

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == '__main__':
    #add all arguments
    parser = ArgumentParser(description='aftermerge v1.0 merging reads based on the mapping information')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-error', help="The error file", type=FileType('w'), dest="errorfile", default=None)
    parser.add_argument('infile', nargs='?', type=FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=FileType('w'), default=sys.stdout)
    #parser.parse_args(['input.txt', 'output.txt']) 
    #Namespace(infile=<_io.TextIOWrapper name='input.txt' encoding='UTF-8'>, outfile=<_io.TextIOWrapper name='output.txt' encoding='UTF-8'>)
    #parser.parse_args([]) 
    #Namespace(infile=<_io.TextIOWrapper name='<stdin>' encoding='UTF-8'>, outfile=<_io.TextIOWrapper name='<stdout>' encoding='UTF-8'>)
    
    #parse the arguments, and set the default output directory 
    args = parser.parse_args()
    #if verbose, print the parameters
    if args.verbose:
        eprint("Verbose: {}".format(args.verbose))
        eprint("Input: {}".format(args.infile))
        eprint("Output: {}".format(args.outfile))
        eprint("Error file: {}".format(args.errorfile))
        
    #infile = open(args.infile.name)
    infile = args.infile
    outfile = args.outfile
    if (args.errorfile is not None):
        errorfile = args.errorfile
    #outfile = open(args.outfile.name, 'w')
    unmapped_reads = 0
    mapped_reads = 0
    merged_pairs = 0
    total_reads = 0
    previousLines = {}
    previousSamReads = {}
    different=0
    processed_lines=0
    for line in infile:
        processed_lines=processed_lines+1
        if (args.verbose and processed_lines%100000==0):
            eprint("Processed {} lines".format(processed_lines))
        line = line.rstrip()
        outline = None
        if (line.startswith("@")):
            #header
            outline = line
        else:
            total_reads = total_reads + 1
            #is read
            samread = samRead.samRead(line)
            #counting for the statistics
            if (samread.ismapped()):
                mapped_reads = mapped_reads + 1
            else:
                unmapped_reads = unmapped_reads + 1
            if (not samread.ispair()):
                #not a paired read
                outline = samread.simpleline
            else:
                #paired read
                if (not samread.pairmapped()):
                    #pair not mapped
                    outline = samread.simpleline
                elif (samread.isSecondaryAlignment()):
                    #read is secondary
                    outline = samread.simpleline
                else:
                    #read can be overlapping
                    if (not samread.qname in previousSamReads):
                        #mate of this read not seen yet
                        previousSamReads[samread.qname] = samread
                    else:
                        #mate seen, so start processing
                        materead = previousSamReads[samread.qname]
                        materead.__class__ = samRead.samRead
                        #remove from dict, to spare memory
                        del previousSamReads[samread.qname]
                        #start checking if overlapping
                        if (not samread.isoverlapping(materead)):
                            #no overlap, so do not fix
                            outline = samread.simpleline + "\n" + materead.simpleline
                        else:
                            #overlapping, so fix
                            #check which read is the first on the reference
                            if (samread.pos > materead.pos):
                                change = samread
                                samread = materead
                                materead = change
                            #check if the cigar string of the overlap is the same,
                            #and get the start of the overlap on the sequence, reference, cigar, matecigar, and length of the overlap
                            (sameoverlapcigar, start, cigarstart, matecigarstart, length, refstart) = samread.overlapHasSameCigar(materead)
                            if (sameoverlapcigar):
                                #same cigar
                                merged_pairs = merged_pairs + 1
                                #get the overlapping sequence, quality and cigar for the first read
                                seq1 = samread.seq[start:start+length]
                                qual1 = samread.qual[start:start+length]
                                mcigar1 = samread.longcigar()[cigarstart:cigarstart+length]
                                #get the overlapping sequence, quality and cigar for the second read
                                start2 = materead.startOfSeqOnRef()
                                seq2 = materead.seq[start2:start2+length]
                                qual2 = materead.qual[start2:start2+length]
                                mcigar2 = materead.longcigar()[matecigarstart:matecigarstart+length]

                                #get the new sequence for the overlap
                                (newseq, newqual, newmcigar) = mergeSequences(seq1, seq2, qual1, qual2, mcigar1, mcigar2, refstart, materead.pos)
                                
                                #correct the new sequence, by adding the non overlapping parts
                                newseq = "{}{}{}".format(samread.seq[0:start], newseq, materead.seq[start2+length:])
                                newqual = "{}{}{}".format(samread.qual[0:start], newqual, materead.qual[start2+length:])
                                newmcigar = "{}{}{}".format(samread.longcigar()[0:cigarstart], newmcigar, materead.longcigar()[matecigarstart+length:])
                                #set flag on 0: forward mapped
                                flag = 0
                                #mapping quality is the mean of both mappings
                                mapq = int((samread.mapq + materead.mapq)/2)
                                cigar = samread.shortcigar(newmcigar)
                                outline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(samread.qname, flag, samread.rname, samread.pos, mapq, cigar, "*", "0", "0", newseq, newqual)
                                
                            else:
                                #different cigar, so overlap is probably an error => add line in error file
                                errorline = "{}\t{}\t{}\t{}".format(samread.rname, refstart, refstart+length, samread.qname)
                                if (args.errorfile is not None):
                                    errorfile.write(errorline + "\n")
                                different=different+1
                                outline = samread.line + "\n" + materead.line
        
        if (outline is not None):
            outfile.write(outline + "\n")
            
    for leftOverReadName in previousSamReads:
        leftOverRead = previousSamReads[leftOverReadName]
        leftOverRead.__class__ = samRead.samRead
        outfile.write(leftOverRead.simpleline + "\n")
            
    infile.close()
    outfile.close()
    if (args.errorfile is not None):
        errorfile.close()
    
        
    #print the statistics to the error output
    eprint("Statistics:")
    eprint("{}\t{} ({}%)".format("Total reads", total_reads, (total_reads/total_reads)*100))
    eprint("{}\t{} ({}%)".format("Unmapped reads", unmapped_reads, (unmapped_reads/total_reads)*100))
    eprint("{}\t{} ({}%)".format("Mapped reads", mapped_reads, (mapped_reads/total_reads)*100))
    eprint("{}\t{} ({}%)".format("Merged pairs", merged_pairs, ((merged_pairs*2)/total_reads)*100))
    eprint("{}\t{} ({}%)".format("Failed merges pairs", different, ((different*2)/total_reads)*100))
    
     
    
    
    
    
