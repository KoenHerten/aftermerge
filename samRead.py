#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 09:51:34 2017

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


class samRead:
    
    def __init__(self, line):
        '''
        initiate all variables
        '''
        self._line = line
        self._linearray = line.split("\t")
        self._binflag = self._tobin()
        
    @property
    def line(self):
        return str(self._line)
    
    @property
    def qname(self):
        return str(self._linearray[0])
        
    @property
    def flag(self):
        return int(self._linearray[1])
        
    @property
    def rname(self):
        return str(self._linearray[2])
        
    @property
    def pos(self):
        return int(self._linearray[3])
        
    @property
    def mapq(self):
        return int(self._linearray[4])
        
    @property
    def cigar(self):
        return str(self._linearray[5])
        
    @property
    def rnext(self):
        return str(self._linearray[6])
        
    @property
    def pnext(self):
        return str(self._linearray[7])
        
    @property
    def tlen(self):
        return int(self._linearray[8])
    
    @property
    def seq(self):
        return str(self._linearray[9])
        
    @property
    def qual(self):
        return str(self._linearray[10])
        
    @property
    def editDistance(self):
        editdist = -1
        if (len(self._linearray) <= 11):
            return -1
        for i in range(11, len(self._linearray)):
            field = str(self._linearray[i])
            if (field.startswith("NM:")):
                fieldArray = field.split(":")
                editdist = fieldArray[2]
        return int(editdist)
        
    @property
    def mismatchPositions(self):
        mispos = None
        if (len(self._linearray) <= 11):
            return None
        for i in range(11, len(self._linearray)):
            field = str(self._linearray[i])
            if (field.startswith("MD:")):
                fieldArray = field.split(":")
                mispos = str(fieldArray[2])
        return mispos
        
    @property
    def simpleline(self):
        simpleline = "{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}".format(
            name=self.qname, flag=self.flag, rname=self.rname, pos=self.pos, mapq=self.mapq, cigar=self.cigar, 
            rnext=self.rnext, pnext=self.pnext, tlen=self.tlen, seq=self.seq, qual=self.qual)
        return str(simpleline)
            
        
    def _tobin(self):
        '''
        function to change the flag to a binary coding
        '''
        flag = int(self.flag)
        r = ""
        while(int(flag) is not 0):
            if (flag % 2 == 0):
                r = "{}{}".format("0", r)
            else:
                r = "{}{}".format("1",r)
            flag = int(flag / 2)
        i=len(r)
        while(i<9):
            r= "0{}".format(r)
            i = i+1
        return r
            
    def ispair(self):
        '''
        function return True if read is part of a pair
        '''
        return (self.isfirst() or self.issecond())
        
    def pairmapped(self):
        '''
        function return True if both reads are mapped
        '''
        return (self.ismapped() and self.ismatemapped())
        
    def isSecondaryAlignment(self):
        '''
        function return True if read is secondary alignment
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-9]) == 1)
            
    def hasMultipleSegments(self):
        '''
        function return True if the read has multiple segments
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-1]) == 1)
            
    
    def ismapped(self):
        '''
        function return True if this read is mapped
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-3]) == 0)
    
    
    def ismatemapped(self):
        '''
        function returns True if pair read is mapped
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-4]) == 0)
    
    def isfirst(self):
        '''
        function returns True if this is the first read in the pair
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-7]) == 1)
            
    def issecond(self):
        '''
        function returns True if this is the second read in the pair
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-8]) == 1)
    
    def isreverse(self):
        '''
        function returns True if this read is reverse on the reference genome
        '''
        binflag = str(self._binflag)
        return (int(binflag[len(binflag)-5]) == 1)
    
    def longcigar(self):
        '''
        function to change the cigar string to a long version ex.: 5M becomes MMMMM
        '''
        cigar = str(self.cigar)
        lcigar = ""
        b = ""
        c=0
        while (c<len(cigar)):
            a = str(cigar[c])
            #if(char ~ /^[0-9]$/):
            if(a.isdigit()):
                b = "{}{}".format(b, a)
            else:
                i=0
                while(i<int(b)):
                    lcigar = "{}{}".format(lcigar, a)
                    i = i + 1
                b=""
            c=c+1
        return lcigar
        
        
    def _mismatchString(self):
        '''
        returns the mismatch positions as long string,
        M for matches, the base in uppercase for mismatches, in lower case for deletions
        '''
        md = str(self.mismatchPositions())
        mstring = ""
        b = ""
        c=0
        isdel = False
        while (c<len(md)):
            a = str(md[c])
            #if(char ~ /^[0-9]$/):
            if(a.isdigit()):
                isdel = False
                b = "{}{}".format(b, a)
            else:
                i=0
                while(i<int(b)):
                    mstring = mstring + "M"
                    i = i + 1
                b=""
                if (a == "^"):
                    isdel = True
                else:
                    base = a.upper()
                    if (isdel):
                        base = a.lower()
                    mstring = mstring + base
            c=c+1
        if (not b == ""):
            i=0
            while(i<int(b)):
                mstring = mstring + "M"
                i = i + 1
            b=""
        return mstring
        
        
    def mismatchcigar(self):
        '''
        returns a long cigar like string,
            soft clipping is s
            match is M
            mismatch is base in upper
            deletion is base in lower
            insertion is i
        '''
        mstring = self._mismatchString()
        lcigar = self.longcigar()
        mcigar = ""
        for c in lcigar:
            if (c == "S"):
                mcigar = mcigar + "s"
            elif (c == "H"):
                #hard clipping: do nothing
                mcigar = mcigar
            else:
                if (c == "I"):
                    mcigar = mcigar + "i"
                else:
                    #match, mismatch or deletion => are in mstring
                    mcigar = mcigar + mstring[0]
                    mstring = mstring[1:]
        return mcigar
    
        
    def shortcigar(self, cigar):
        '''
        function to change the long cigar string to a normal one ex.: MMMMM becomes 5M
        '''
        cigar = str(cigar)
        scigar = ""
        c = ""
        l = ""
        for letter in cigar:
            if (l == letter):
                c = c + 1
            else:
                scigar = "{}{}{}".format(scigar, c, l)
                c = 1
                l = letter
        scigar = "{}{}{}".format(scigar, c, l)
        return scigar
        
    def getLengthOnReference(self):
        '''
        returns the length the sequence on the reference
        (counting all M and Ds in the cigar string)
        '''
        lcigar = str(self.longcigar())
        return lcigar.count('M') + lcigar.count('D')
        
    def isoverlapping(self, samread):
        '''
        function return True if both reads are overlapping
        '''
        if (not isinstance(samread, samRead)):
            return False
        if (not self.rname == samread.rname):
            #different reference
            return False
        else:
            #same reference
            if(self.pos <= samread.pos):
                #this is before samread
                if (self.pos + self.getLengthOnReference() -1 >= samread.pos):
                    #samread starts in the mapped area of this
                    return True
                else:
                    #samreads start out of this area
                    return False
            else:
                #samread is before this
                if (samread.pos + samread.getLengthOnReference() -1 >= self.pos):
                    #this starts in the mapping area of samread
                    return True
                else:
                    #this start out of samreads area
                    return  False
        
    def overlapHasSameCigar(self, samread):
        '''
        function returns (true/false, start)
            True if the cigar of the overlap is matching
            start the start in the first sequence
        '''
        if (not isinstance(samread, samRead)):
            return (False, -1, None, None, None, None)
        if (not self.isoverlapping(samread)):
            return (False, -1, None, None, None, None)
        #overlapping reads
        firstread = self
        secondread = samread
        if (secondread.pos > firstread.pos):
            change = secondread
            secondread = firstread
            firstread = change
        seqstart = firstread.pos - secondread.pos
        cigarstart = seqstart
        samlcigar = str(secondread.longcigar())
        samlcigar = samlcigar.strip("S")
        #TODO include I and D into calculations
        changestart = True
        insertions = 0
        while(changestart):
            newinsertions = samlcigar[0:cigarstart].count('I')
            mcount = samlcigar[0:cigarstart].count('M')
            if (not mcount == (firstread.pos - secondread.pos)):
                seqstart = seqstart + 1
                cigarstart = cigarstart + 1
            if(newinsertions == insertions):
                #same number of insertions, so check if next base is not an insertion
                if (len(samlcigar) > seqstart):
                    if (samlcigar[seqstart] == "M"):
                        changestart = False
                    elif (samlcigar[seqstart] == "I"):
                        newinsertions = newinsertions + 1
                        seqstart = seqstart + 1
                        cigarstart = cigarstart + 1
                    elif (samlcigar[seqstart] == "D"):
                        #have a delition
                        cigarstart = cigarstart + 1
                        changestart = False
                else:
                    changestart = False
            else:
                #different number of insertions
                insertions = newinsertions
                seqstart = seqstart
                cigarstart = cigarstart
            
        isstart=True
        for letter in str(secondread.longcigar()):
            if (isstart and letter is "S"):
                cigarstart = cigarstart + 1
            else:
                isstart = False
        length = min((secondread.getLengthOnReference() - cigarstart), firstread.getLengthOnReference())
        samlcigar = str(secondread.longcigar())[cigarstart:cigarstart+length]
        matelcigar = str(firstread.longcigar())
        matecigarstart = 0
        isstart=True
        for letter in matelcigar:
            if (isstart and letter is "S"):
                matecigarstart = matecigarstart + 1
            else:
                isstart = False
        matelcigar = matelcigar[matecigarstart:matecigarstart+length]
        refstart = secondread.pos + str(secondread.longcigar())[0:cigarstart].count('M') + str(secondread.longcigar())[0:cigarstart].count('D')
        if (str(samlcigar) == str(matelcigar)):
            return (True, seqstart, cigarstart, matecigarstart, length, refstart)
        else:
            if (length == -1):
                #not overlapping, but has to be connected
                return (True, seqstart, cigarstart, matecigarstart, length, refstart)
            return (False, seqstart, cigarstart, matecigarstart, length, refstart)
            
    def startOfSeqOnRef(self):
        '''
        returns the index of the sequence, which is the first base on the reference
        '''
        i = 0
        longcigar = self.longcigar()
        while (i < len(longcigar)):
            if (longcigar[i] != "S"):
                return i
            i = i +1
        return i
        
    def endOfSeqOnRef(self):
        '''
        returns the index of the sequence, which is the last base on the reference
        '''
        longcigar = self.longcigar()
        i = len(longcigar) -1
        while (i > 0):
            if (longcigar[i] != "S"):
                return i
            i = i -1
        return i
        
        
