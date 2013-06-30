#!/usr/bin/env python

import sys
import csv
import subprocess
from collections import namedtuple
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def overlap(pos1,pos2):
    (start1,end1)=map(int,pos1)
    (start2,end2)=map(int,pos2)
    orange='-'
    if ((start1 <= end2) and (start2 <= end1)):   
        max_start = max(start1, start2)
        min_end = min(end1, end2)
        orange='-'.join([str(min(start1, start2)),str(max(end1, end2))])
    return orange

def readblast(file):
    #Read blast results
    recblast = namedtuple('recblast','qw, sub, pid, len, mism, gaps, qstart, qend, sstart, send, eval, score')
   # nbhits=0
    blast=[]
    for line in csv.reader(open(sys.argv[1],'rb'),delimiter='\t'):
        record=recblast._make(line)
        #print record
        blast.append(record)
    return blast

def overlapping(list):
    coord=[]
    for pos in list:
        coord.append(('S',pos[0]))
        coord.append(('E',pos[1]))
    coord.sort(key = lambda x : x[0], reverse = True)
    coord.sort(key = lambda x : x[1])

    spos=0
    ct=0
    ovl=[]
    #print coord
    for pos in coord:
        if pos[0]=='S':
            ct+=1
            if ct==2:
                spos=pos[1]
        if pos[0]=='E':
            ct-=1
            if ct==0:
                ovl.append((spos,pos[1]))
    return ovl
    

def retr_db(database,sid):
    child = subprocess.Popen("blastdbcmd -db %s -entry '%s'" % (database,sid),stdout=subprocess.PIPE,shell=True)
    (out,err)=child.communicate()
    allseq=''.join(out.split('\n')[1:])
    return allseq



###MAIN

    
blast=readblast(sys.argv[1])
Hits=defaultdict(list)
for rec in blast:
    Hits[rec.sub].append(rec)
    #qgene=rec.qw.split('_')[1]

print len(blast),"blast hits in",len(Hits),"scaffold"

seHits=[]

for loc in Hits:
    posn=[]
    for rec in Hits[loc]:
        spos=sorted([int(rec.sstart),int(rec.send)])
        posn.append((spos[0],spos[1])) 
        
    #print posn
    #recover not overlapping hit positions on scaffold
    redposn=overlapping(posn)
    ##gather all records associated to a hit position
    #print loc,len(posn),'hits',len(redposn),'hit positions:'
    loc_hits=defaultdict(list)
    for pos in redposn:
        for rec in Hits[loc]:
            if not overlap((rec.sstart,rec.send),pos)=='-':
                loc_hits[pos].append(rec)
    # select best hit for each hit position
    bscore=float(0)
    for pos in loc_hits:
        brec=loc_hits[pos][0]
        for rec in loc_hits[pos]:
            #print pos,rec.qw,rec.score
            if float(rec.score)>float(brec.score):
                brec=rec
        seHits.append((loc,pos,brec.qw,brec))

seHits.sort(key = lambda x : x[2])

qlen=60
database=sys.argv[2]
prefix=sys.argv[3]
prot=open(sys.argv[1].split('.')[0]+'.prot.fa','w')
nucl=open(sys.argv[1].split('.')[0]+'.nucl.fa','w')
scaf=open(sys.argv[1].split('.')[0]+'.scaf.fa','w')
for elt in seHits:
    rec=elt[3]
    sens='+' if rec.sstart<rec.send else '-'
    (sstart,send)=map(int,sorted([rec.sstart,rec.send],key=int))
    (qstart,qend)=map(int,sorted([rec.qstart,rec.qend],key=int))
    rawseq=retr_db(database, elt[0])
    clen=len(rawseq)
    scaf.write('>'+elt[0]+'_'+rec[2]+'\n'+rawseq+'\n')
    if sens=='+':
        start=sstart-3*(qstart-1) if sstart-3*(qstart-1)>1 else sstart
        end=send+3*(qlen-qend) if send+3*(qlen-qend)<clen else send  
    elif sens=='-':
        start=sstart-3*(qlen-qend) if sstart-3*(qlen-qend)>1 else sstart
        end=send+3*(qstart-1) if send+3*(qstart-1)<clen else send  
    
    print elt[0],elt[2],sens,sstart,send,qstart,qend,start,end,clen
    
    hdseq=Seq(rawseq[start-1:end],IUPAC.unambiguous_dna)
    if sens=='-':
        hdseq=hdseq.reverse_complement()
    hdseq=hdseq.split('N')[0]
    print hdseq.translate()
    #head='>'+prefix+'_'+elt[2].split('_')[1]+'-'+elt[0].split('|')[1]+'@'+'-'.join(map(str,elt[1]))+'\n'    
    head='>'+prefix+'_'+elt[2].split('_')[1]+'|ND'+elt[0].split('_')[1]+'@'+'-'.join(map(str,elt[1]))+'\n'
    prot.write(head+hdseq.translate().tostring()+'\n')
    nucl.write(head+hdseq.tostring()+'\n')

#+'-'.join(map(str,elt[1])







