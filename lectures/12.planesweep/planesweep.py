#!/usr/local/canopy/bin/python

## sudo canopy_cli setup /usr/local/canopy/

PLOT_FIGURES = 1


## Plane-sweep coverage analysis
###############################################################################
##
## In this notebook, we will consider the problem of computing the coverage
## across a sequence given a set of alignments. The alignments could
## come from reads aligned to a genome, or the layout of reads within a newly
## assembled contig. It could also be used to track arbitrary features, 
## such as the gene regions or other annotations.
## 
## For example, given these five read alignments, the goal is to compute the
## coverage vector below it
##
##
## pos:    1    10    20   30   40   50    60    70    80     90
## r1:     [==========================]
## r2:          [================]
## r3:                 [================================]
## r4:                      [=====================]
## r5:                                     [===================]
##
## cov:    11111222222233333444444333332222333333332222221111111
##
## Note the algorithms and data structures presented are also generally
## useful for many other genomic analyses working with intervals
##
###############################################################################






## Load layout information
###############################################################################
## 
## To get started, we will first load in the layout of PacBio reads within 
## a newly assembled contig built using the Celera Assemler. The information 
## of how reads are placed is stored within the "tigStore", which is an on-disk 
## binary database. Because it is a binary database, you cannot read it directly,
## but instead should use the 'tigStore' command to query or update the database.
##
## These commands will display the unitig layouts:
##
## $ cd /seq/schatz/mnattest/skbr3/local_assembly/assembly_HLA_Feb17/CA/r1
## 
## $ /sonas-hs/schatz/hpc/home/gurtowsk/sources/wgs.svn/Linux-amd64/bin/tigStore -g r1.gkpStore -t r1.tigStore/ 1 -d layout -U | head -20
## 
## unitig 0
## len 0
## cns
## qlt
## data.unitig_coverage_stat  1.000000
## data.unitig_microhet_prob  1.000000
## data.unitig_status         X
## data.unitig_suggest_repeat F
## data.unitig_suggest_unique F
## data.unitig_force_repeat   F
## data.unitig_force_unique   F
## data.contig_status         U
## data.num_frags             2884
## data.num_unitigs           0
## FRG type R ident     18318 container         0 parent     16010 hang   2962   1054 position  19674      0
## FRG type R ident      2248 container     18318 parent     18318 hang    217  -4192 position  15480    218
## FRG type R ident      8536 container      2248 parent      2248 hang    609  -4597 position    827  10881
## FRG type R ident     10602 container      2248 parent      2248 hang    848   -420 position  15057   1067
## FRG type R ident     16010 container         0 parent     18318 hang   1054   2962 position  22636   1089
## FRG type R ident      4611 container     18318 parent     18318 hang   1155  -3014 position   1157  16657
##
## This shows the first few reads (fragments, coded as FRG) within unitig 0.
## The most important values are the last two numbers for the FRG records which 
## encodes the position of the read in the contig.
## 
## For example, the first read in this unitig is read 18318, and spans from
## position 19674 to position 0. Because the start position is greater than
## the end, we know this read has been reverse complemented and really starts
## at position 0.
##
## For PacBio assemblies, we are often interested in how well the longest reads
## were assembled. This command will determine which unitig has the longest read in it:
## 
## $ /sonas-hs/schatz/hpc/home/gurtowsk/sources/wgs.svn/Linux-amd64/bin/tigStore -g r1.gkpStore -t r1.tigStore/ 1 -d layout -U | awk '{if (/^unitig/){tig=$2}else if(/^FRG/){l=$15-$14;if(l<0){l*=-1;}
## 
## 3   37839  FRG  type  R  ident  5708   container  0  parent  10683  hang  122   20199  position  416262  454101
## 3   36519  FRG  type  R  ident  5022   container  0  parent  167    hang  1199  9553   position  402020  365501
## 4   35962  FRG  type  R  ident  3652   container  0  parent  10438  hang  199   13564  position  15094   51056
## 0   34550  FRG  type  R  ident  11392  container  0  parent  11328  hang  1010  18989  position  182086  216636
## 15  34302  FRG  type  R  ident  2622   container  0  parent  5938   hang  6855  7174   position  948470  914168
## 15  33695  FRG  type  R  ident  16040  container  0  parent  16226  hang  5741  9003   position  539963  573658
## 15  33671  FRG  type  R  ident  5938   container  0  parent  11492  hang  3720  16069  position  941296  907625
## 17  33504  FRG  type  R  ident  2341   container  0  parent  9030   hang  192   15632  position  382387  415891
## 3   33148  FRG  type  R  ident  8132   container  0  parent  7619   hang  1171  14849  position  162425  129277
## 3   32963  FRG  type  R  ident  5377   container  0  parent  14108  hang  1184  17010  position  102322  69359
## 
## This means that unitig 3 (column 1) has the longest read at 38,839bp which is read 5708. 
## It also has the second longest read with 36,519bp, which is read 5022.
##
## This command will report the read layout for just unitig 3, and just save the 
## readid, start position, and end position to a file:
## 
## $ /sonas-hs/schatz/hpc/home/gurtowsk/sources/wgs.svn/Linux-amd64/bin/tigStore -g r1.gkpStore -t r1.tigStore/ 1 -d layout -u 3 | grep '^FRG' | awk '{print $5,$14,$15}' > ~/readid.start.stop.txt
## 
## $ head -3 ~/readid.start.stop.txt
## 1 0 19814
## 2 799 19947
## 3 1844 13454
## 
## $ tail -3 ~/readid.start.stop.txt
## 1871 973590 965902
## 1872 966703 973521
## 1873 973632 966946
## 
## Now we can load this file and examine the coverage across the contig
##
###############################################################################





## Load the read placement file
###############################################################################
##
## File format is: readid \t startpos \t endpos
##                 If startpos > endpos, the read has been reverse complemented
## 
## It would also be very easy to generate this type of file from a SAM or BED file
##
###############################################################################


import matplotlib.pyplot as plt
from collections import deque
import heapq
import time


## Limit the number of reads to load, -1 for unlimited
MAXREADS = -1

## Limit the number of reads to plot
MAX_READS_LAYOUT = 500

## Path to reads file
READFILE = "readid.start.stop.txt"
print "Loading reads from " + READFILE
f = open(READFILE)

## This will have a list of tuples containing (readid, start, end, rc)
## rc is 0 if the read was originally oriented forward, 1 for reverse

reads = []
totallen = 0
readidx = 0

for line in f:
  line = line.rstrip()
  fields = line.split()
  readid  = fields[0]
  start   = int(fields[1])
  end     = int(fields[2])
  rc      = 0

  if (start > end):
    start = int(fields[2])
    end   = int(fields[1])
    rc = 1

  readinfo = (readid, start, end, rc)
  reads.append(readinfo)

  if (end > totallen):
    totallen = end

  readidx += 1
  if (readidx == MAXREADS):
    break

print "Loaded layout information for %d reads" % (len(reads))

for i in xrange(3):
  print "  %d: %s [%d - %d] %d" % (i, reads[i][0], reads[i][1], reads[i][2], reads[i][3])

print "  ..."

for i in xrange(len(reads)-3, len(reads)):
  print "  %d: %s [%d - %d] %d" % (i, reads[i][0], reads[i][1], reads[i][2], reads[i][3])

print "\n\n"



## Plot the read layouts
###############################################################################

PLOT_LAYOUT = 0 + PLOT_FIGURES

if (PLOT_LAYOUT):

  plt.figure()

  print "Plotting layout\n\n"

  ## draw the layout of reads
  for i in xrange(min(MAX_READS_LAYOUT, len(reads))):
    r = reads[i]
    readid = r[0]
    start  = r[1]
    end    = r[2]
    rc     = r[3]
    color  = "blue"

    if (rc == 1): 
      color = "red"

    plt.plot ([start,end], [-2*i, -2*i], lw=4, color=color)

  plt.draw()
  plt.show()



## Brute force the coverage profile
###############################################################################
##
## There are several techniques for computing the coverage profile. One
## basic approach would be to separately compute how many reads span each 
## position in the genome. This would take O(G * N) time where G is the 
## length of the genome, and N is the number of reads since it requires O(1)
## time to detemine if a read spans a given position:
##
## if ((startpos <= querypos) and (querypos <= endpos)):
##    print "This read spans the query position!"
##
## Here we do something slightly smarter which is to initialize the coverage
## array with all zeros, and then process each read to increment the coverage
## array at each position in the read. This requires O(G + N*L) where L 
## is the average length of the read. Since L << G, this can be a huge advantage
##
###############################################################################

print "Brute force computing coverage over %d bp" % (totallen)

starttime = time.time()
brutecov = [0] * totallen

for r in reads:
 # print " -- [%d, %d]" % (r[1], r[2])

  for i in xrange(r[1], r[2]):
    brutecov[i] += 1

brutetime = (time.time() - starttime) * 1000.0

print "  Brute force complete in %0.02f ms" % (brutetime)
print brutecov[0:10]


## Compute the max coverage
###############################################################################
##
## With the coverage profile, it is trivial to compute the max coverage or
## other simple queries
##
###############################################################################

maxcov = 0
lowcov = 0
lowcovthresh = 10

for c in brutecov:
  if c > maxcov:
    maxcov = c
  if c <= lowcovthresh:
    lowcov += 1

print "max cov: %d, there are %d low coverage bases (<= %d depth)" % (maxcov, lowcov, lowcovthresh)

print "\n\n"


## Show the read depth histogram
###############################################################################

PLOT_DEPTH_HISTOGRAM = 0 + PLOT_FIGURES

if (PLOT_DEPTH_HISTOGRAM):
  plt.figure()
  plt.hist(brutecov)
  plt.show()



## Delta encode the coverage profile
###############################################################################
##
## Imagine the coverage profile looks like this:
##
##         ##                                    #####         ###
##        ####                               ############### #########
##    ############  #########  ###       ################### ########### 
##  ###############################   ###################### #############
##
##  11222234432222112222222221122210001112222333344444333333033444333322110
##  
##  00000000001111111111222222222233333333334444444444555555555566666666667
##  01234567890123456789012345678901234567890123456789012345678901234567890
##
## 
##  Most consecutive positions have the same depth 
##     --> we can recover exactly this plot if we just record the transitions
##
##         #                                     #             #  
##        ## #                               #   #    #      # #  #   
##    #   ## ##     #          #         #   #   #    #      # #  #   #  
##  # #   ## ##   # #        # #  #   #  #   #   #    #      # #  #   # # 
##
##  1 2   34 32   1 2        1 2  10  1  2   3   4    3      3 4  3   2 1 0 
##
##  0 0   00 01   1 1        2 2  33  3  3   4   4    5     55 5  6   6 6 7
##  0 2   67 90   4 6        5 7  01  4  7   1   5    0     67 9  2   6 8 0
##
##  Instead of 70 depth values, we can record just 24 without loss of information!
## 
###############################################################################

print "Delta encoding coverage plot"

starttime = time.time()

deltacov = []
curcov = -1
for i in xrange(0, len(brutecov)):
  if brutecov[i] != curcov:
    curcov = brutecov[i]
    delta = (i, curcov)
    deltacov.append(delta)

## Finish up with the last position
deltacov.append((totallen, 0))

print "Delta encoding required only %d steps, saving %0.02f%% of the space in %0.02f ms" % (len(deltacov), (100.0*float(totallen-len(deltacov))/totallen), (time.time()-starttime) * 1000.0)

for i in xrange(3):
  print "  %d: [%d,%d]" % (i, deltacov[i][0], deltacov[i][1])

print "  ..."

for i in xrange(len(deltacov)-3, len(deltacov)):
  print "  %d: [%d,%d]" % (i, deltacov[i][0], deltacov[i][1])

print "\n\n"

## Plot the coverage profile
###############################################################################

PLOT_COVERAGE_PROFILE = 0 + PLOT_FIGURES

if (PLOT_COVERAGE_PROFILE):

  plt.figure()

  print "Plotting coverage profile\n\n"

  ## draw the layout of reads
  for i in xrange(min(MAX_READS_LAYOUT, len(reads))):
    r = reads[i]
    readid = r[0]
    start  = r[1]
    end    = r[2]
    rc     = r[3]
    color  = "blue"

    if (rc == 1): 
      color = "red"

    plt.plot ([start,end], [-2*i, -2*i], lw=4, color=color)

  ## draw the base of the coverage plot
  plt.plot([0, totallen], [0,0], color="black")

  ## draw the coverage plot
  for i in xrange(len(deltacov)-1):
    x1 = deltacov[i][0]
    x2 = deltacov[i+1][0]
    y1  = deltacov[i][1]
    y2  = deltacov[i+1][1]

    ## draw the horizonal line
    plt.plot([x1, x2], [y1, y1], color="black")
    
    ## and now the right vertical to the new coverage level
    plt.plot([x2, x2], [y1, y2], color="black")

  plt.draw()
  plt.show()



## Plane-Sweep Coverage
###############################################################################
##
## Notice the coverage only changes at the beginning or end of a read
## so just need to walk from read to read and keep track of depth along the way
##
## Returning to our example again:
##
##  pos:    1    10    20   30   40   50    60    70    80     90
##  r1:     [==========================]
##  r2:          [================]
##  r3:                 [================================]
##  r4:                      [=====================]
##  r5:                                     [===================]
##
##  cov:    11111222222233333444444333332222333333332222221111111
##
## We can imagine scanning across this from left to right to
## keep track of how many reads span each position. This is a widely
## used approach in computational geometry called a plane-sweep
## algorithm.
##
## The basic technique follows like this:
##   - Assume layout is in sorted order by start position
##   - walking from start position to start position
##       - check to see if we past any read ends 
##            - coverage goes down by one when a read ends
##       - coverage goes up by one when new read is encountered
##
## Core concept: 
##   - use a list/deque/heap to track how many reads currently intersect the plane
##   - the number of elements in the list corresponds to the current depth
##
## In this example, the plane sweeps like this:
##
##     1 (add 50): 50
##    10 (add 40): 40, 50 <- notice insert out of order 
##    20 (add 80): 40, 50, 80
##    30 (add 70): 40, 50, 70, 80 <- out of order again
##    40 (sub 40): 50, 70, 80
##    50 (sub 50): 70, 80
##    60 (add 90): 70, 80, 90
##    70 (sub 70): 80, 90
##    80 (sub 80): 90
##    90 (sub 90): 0
##
###############################################################################


print "Beginning list-based plane sweep over %d reads" % (len(reads))

starttime = time.time()

## record the delta encoded depth using a plane sweep
deltacovplane = []

## use a list to record the end positions of the elements currently in plane
planelist = []

## BEGIN SWEEP 
for r in reads:
  startpos = r[1]
  endpos   = r[2]

  ## clear out any positions from the plane that we have already moved past
  while (len(planelist) > 0):

    if (planelist[0] <= startpos):
      ## the coverage steps down, extract it from the front of the list 
      oldend = planelist.pop(0) 
      deltacovplane.append((oldend, len(planelist)))
    else:
      break

  ## Now insert the current endpos into the correct position into the list
  insertpos = -1
  for i in xrange(len(planelist)):
    if (endpos < planelist[i]):
        insertpos = i
        break

  if (insertpos > 0):
    planelist.insert(insertpos, endpos)
  else:
    planelist.append(endpos)


  ## Finally record that the coverage has increased
  deltacovplane.append((startpos, len(planelist)))


## Flush any remaining end positions
while (len(planelist) > 0):
  oldend = planelist.pop(0) 
  deltacovplane.append((oldend, len(planelist)))


## Report statistics
planelisttime = (time.time() - starttime) * 1000.0
print "Plane sweep found %d steps, saving %0.02f%% of the space in %0.2f ms (%0.02f speedup)!" % (len(deltacovplane), (100.0*float(totallen-len(deltacovplane))/totallen), planelisttime, brutetime/planelisttime)

for i in xrange(3):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "  ..."

for i in xrange(len(deltacovplane)-3, len(deltacovplane)):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "\n\n"



## Checking for errors
###############################################################################
##
## The above implementation is flawed because there may be multiple events 
## (read starts and/or read ends) at the same position. This will leading to 
## multiple entries at the same position, although there should only be one
##
###############################################################################

print "Checking for duplicates:"

lastpos = -1
for i in xrange(len(deltacovplane)):
  if deltacovplane[i][0] == lastpos:
    print "  Found duplicate: %d: [%d,%d] == %d: [%d,%d]" % (i-1, deltacovplane[i-1][0], deltacovplane[i-1][1], i, deltacovplane[i][0], deltacovplane[i][1])
  lastpos = deltacovplane[i][0]

print "\n\n"





## Fix the update rule
###############################################################################
##
## Lets fix the implementation to check for multiple events at the same position
## by peeking ahead to check for multiple ends or multiple starts at the exact
## same position
## 
## Strictly speaking this may give slightly different results than the brute
## force coverage profile since this will create an event if a read ends
## at exactly the same place another begins, whereas the brute force approach
## will link together those regions because they have the same coverage
##
###############################################################################

print "Beginning correct list-based plane sweep over %d reads" % (len(reads))

starttime = time.time()

## record the delta encoded depth using a plane sweep
deltacovplane = []

## use a list to record the end positions of the elements currently in plane
planelist = []

## BEGIN SWEEP (note change to index based so can peek ahead)
for rr in xrange(len(reads)):
  r = reads[rr]
  startpos = r[1]
  endpos   = r[2]

  ## clear out any positions from the plane that we have already moved past
  while (len(planelist) > 0):

    if (planelist[0] <= startpos):
      ## the coverage steps down, extract it from the front of the list 
      oldend = planelist.pop(0) 

      nextend = -1
      if (len(planelist) > 0):
        nextend = planelist[0]

      ## only record this transition if it is not the same as a start pos
      ## and only if not the same as the next end point
      if ((oldend != startpos) and (oldend != nextend)):
        deltacovplane.append((oldend, len(planelist)))
    else:
      break

  ## Now insert the current endpos into the correct position into the list
  insertpos = -1
  for i in xrange(len(planelist)):
    if (endpos < planelist[i]):
        insertpos = i
        break

  if (insertpos > 0):
    planelist.insert(insertpos, endpos)
  else:
    planelist.append(endpos)


  ## Finally record that the coverage has increased
  ## But make sure the current read does not start at the same position as the next
  if ((rr == len(reads)-1) or (startpos != reads[rr+1][1])):
    deltacovplane.append((startpos, len(planelist)))

  ## if it is at the same place, it will get reported in the next cycle

## Flush any remaining end positions
while (len(planelist) > 0):
  oldend = planelist.pop(0) 
  deltacovplane.append((oldend, len(planelist)))


## Report statistics
planelisttime = (time.time() - starttime) * 1000.0
print "Plane sweep found %d steps, saving %0.02f%% of the space in %0.2f ms (%0.02f speedup)!" % (len(deltacovplane), (100.0*float(totallen-len(deltacovplane))/totallen), planelisttime, brutetime/planelisttime)

for i in xrange(3):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "  ..."

for i in xrange(len(deltacovplane)-3, len(deltacovplane)):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "\n\n"

print "Checking for duplicates in corrected version:"

lastpos = -1
for i in xrange(len(deltacovplane)):
  if deltacovplane[i][0] == lastpos:
    print "  Found duplicate: %d: [%d,%d] == %d: [%d,%d]" % (i-1, deltacovplane[i-1][0], deltacovplane[i-1][1], i, deltacovplane[i][0], deltacovplane[i][1])
  lastpos = deltacovplane[i][0]

print "\n\n"



## Heap-based plane-sweep algorithm
###############################################################################
##
## Now lets replace the list with a better data structure
##
## Lists are bad for maintaining the plane because:
##   removing from the front requires shifting everything over: O(N) time
##   inserting into the middle (or end) also requires shifting and possibly resizing: O(N)
##
## This can result in a worst case runtime of O(N^2)
##  - for each read, shift the list by N positions
##
## Fortunately, with most real datasets the maximum length of the list is bounded by the
## maximum read depth (D), so will result in an O(N*D) runtime
##
## If the reads (intervals) were all the same length, we would always append to the end
## and remove from the very front. A good data structure for this is called a deque
## and enables O(1) remove-from-front and O(1) append-to-end operations
##
## It is implemented as a doublely-linked-list, and can be used from the collections API
##
## from collections import deque
## d = deque()
## d.append()  <- adds to the end
## d.popleft() <- removes from the front
## 
## this estalishes a "queue" for the end positions and accelerates to O(N) runtime!
##
## However, in general, the intervals will not be of uniform length so this will not
## be sufficient. The data structure we really want should support O(1) ability
## to find the minimum item of the plane, "fast" ability to remove it, 
## and "fast" ability to add a new item to the plane
##
## A good data structure for this is called a "heap" (aka min-heap, aka priority queue)
##
## Heaps are:
##  - complete binary trees such that the parent is smaller than both of the children
##     - the height of the tree is completely balanced, except for the bottommost level
##     - the relative ordering of the left and right children is arbitrary
##  - finding the minimum value is trival  O(1)
##     - by construction it is always the top of the heap
##  - removing the minimum requires O(log n) steps: 
##      move a leaf value into the root position
##      if it is smaller than children, done
##      else iteratively swap it with the smaller of the two children
##  - adding a new value is fast: O(log n) worst case, O(1) for sorted values, O(1) average case
##      add the new value into a free cell at the bottom of the heap
##      then up-heap the element until the heap-invariant is met
##
## We can store a heap by embedding the tree inside of an array
##   - the children at node k are indexed at nodes 2k+1 and 2k+2
##   - this way we dont have to maintain extra pointers 
##   - up-heaps and down-heaps are just swapping values in the list
##   - finding the next free leaf is just the first free cell in the array
##
##                          0
##                  1              2
##               3     4         5     6
##              7 8   9 10     11 12 13 14
##
##   =>
## 
##             0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
##
##   The next free cell is 15, which is a child of 7 (2*7+1 = 15)
##
## In python, these techniques are implemented in the heapq package
##
## This restores us to O(N log N) worst case runtime. However, the time to 
## maintain the heap is bounded by the size of it, so in practice the runtime
## will be much faster: O(N log D) where D is the maximum depth of coverage.
##
###############################################################################

print "Beginning heap-based plane sweep over %d reads" % (len(reads))

starttime = time.time()

## record the delta encoded depth using a plane sweep
deltacovplane = []

## use a list to record the end positions of the elements currently in plane
planeheap = []

## BEGIN SWEEP (note change to index based so can peek ahead)
for rr in xrange(len(reads)):
  r = reads[rr]
  startpos = r[1]
  endpos   = r[2]

  ## clear out any positions from the plane that we have already moved past
  while (len(planeheap) > 0):

    if (planeheap[0] <= startpos):
      ## the coverage steps down, extract it from the front of the list 
      ## oldend = planelist.pop(0) 
      oldend = heapq.heappop(planeheap)

      nextend = -1
      if (len(planeheap) > 0):
        nextend = planeheap[0]

      ## only record this transition if it is not the same as a start pos
      ## and only if not the same as the next end point
      if ((oldend != startpos) and (oldend != nextend)):
        deltacovplane.append((oldend, len(planeheap)))
    else:
      break

  ## Now insert the current endpos into the correct position into the list

  ##insertpos = -1
  ##for i in xrange(len(planelist)):
  ##  if (endpos < planelist[i]):
  ##      insertpos = i
  ##      break

  ##if (insertpos > 0):
  ##  planelist.insert(insertpos, endpos)
  ##else:
  ##  planelist.append(endpos)

  heapq.heappush(planeheap, endpos)


  ## Finally record that the coverage has increased
  ## But make sure the current read does not start at the same position as the next
  if ((rr == len(reads)-1) or (startpos != reads[rr+1][1])):
    deltacovplane.append((startpos, len(planeheap)))

  ## if it is at the same place, it will get reported in the next cycle


## Flush any remaining end positions
while (len(planeheap) > 0):
  ##oldend = planelist.pop(0) 
  oldend = heapq.heappop(planeheap)
  deltacovplane.append((oldend, len(planeheap)))



## Report statistics
planeheaptime = (time.time()-starttime) * 1000.0
print "Heap-Plane sweep found %d steps, saving %0.02f%% of the space in %0.02f ms (%0.02f speedup)!" % (len(deltacovplane), (100.0*float(totallen-len(deltacovplane))/totallen), planeheaptime, brutetime/planeheaptime)

for i in xrange(3):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "  ..."

for i in xrange(len(deltacovplane)-3, len(deltacovplane)):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "\n\n"

print "Checking for duplicates in heap-ified version:"

lastpos = -1
for i in xrange(len(deltacovplane)):
  if deltacovplane[i][0] == lastpos:
    print "  Found duplicate: %d: [%d,%d] == %d: [%d,%d]" % (i-1, deltacovplane[i-1][0], deltacovplane[i-1][1], i, deltacovplane[i][0], deltacovplane[i][1])
  lastpos = deltacovplane[i][0]

print "\n\n"




## Max-coverage read ids
###############################################################################
##
## Now lets use this framework to figure out which reads span the 
## position of max coverage. The technique is to update the heap to store
## the read end position as the primary key, but also associate it
## with the read id. Then whenever we see we have established a new maximum
## coverage level, we record the ids of the reads in the heap.
##
###############################################################################



print "Beginning heap-based plane sweep over %d reads to search for max reads" % (len(reads))

starttime = time.time()

## record the delta encoded depth using a plane sweep
deltacovplane = []

## use a list to record the end positions of the elements currently in plane
planeheap = []

## use a list to record the ids at max coverage
maxcov = -1
maxcovpos = -1
maxcovreads = []

## BEGIN SWEEP (note change to index based so can peek ahead)
for rr in xrange(len(reads)):
  r = reads[rr]
  readid   = r[0]
  startpos = r[1]
  endpos   = r[2]

  ## clear out any positions from the plane that we have already moved past
  while (len(planeheap) > 0):

    if (planeheap[0][0] <= startpos):
      ## the coverage steps down, extract it from the front of the list 
      oldend = heapq.heappop(planeheap)[0]

      nextend = -1
      if (len(planeheap) > 0):
        nextend = planeheap[0][0]

      ## only record this transition if it is not the same as a start pos
      ## and only if not the same as the next end point
      if ((oldend != startpos) and (oldend != nextend)):
        deltacovplane.append((oldend, len(planeheap)))
    else:
      break

  ## Now insert the current endpos into the correct position into the list
  heapq.heappush(planeheap, (endpos, readid))


  ## Finally record that the coverage has increased
  ## But make sure the current read does not start at the same position as the next
  if ((rr == len(reads)-1) or (startpos != reads[rr+1][1])):
    cov = len(planeheap)
    if (cov > maxcov):
      maxcov = cov
      maxcovpos = startpos
      maxcovreads = []
      for rr in planeheap:
        maxcovreads.append(rr[1])

    deltacovplane.append((startpos, len(planeheap)))


## Flush any remaining end positions
while (len(planeheap) > 0):
  oldend = heapq.heappop(planeheap)[0]
  deltacovplane.append((oldend, len(planeheap)))



## Report statistics
planeheaptime = (time.time()-starttime) * 1000.0
print "Heap-Plane sweep found %d steps, saving %0.02f%% of the space in %0.02f ms (%0.02f speedup)!" % (len(deltacovplane), (100.0*float(totallen-len(deltacovplane))/totallen), planeheaptime, brutetime/planeheaptime)

print "The %d reads at the position of maximum depth (%d) are:" % (maxcov, maxcovpos)
print maxcovreads

print "\n\n"

for i in xrange(3):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "  ..."

for i in xrange(len(deltacovplane)-3, len(deltacovplane)):
  print "  %d: [%d,%d]" % (i, deltacovplane[i][0], deltacovplane[i][1])

print "\n\n"

print "Checking for duplicates in heap-ified version:"

lastpos = -1
for i in xrange(len(deltacovplane)):
  if deltacovplane[i][0] == lastpos:
    print "  Found duplicate: %d: [%d,%d] == %d: [%d,%d]" % (i-1, deltacovplane[i-1][0], deltacovplane[i-1][1], i, deltacovplane[i][0], deltacovplane[i][1])
  lastpos = deltacovplane[i][0]

print "\n\n"




## Plot the max depth position
###############################################################################

PLOT_MAX_DEPTH = 0 + PLOT_FIGURES

if (PLOT_MAX_DEPTH):

  plt.figure()

  print "Plotting max depth\n\n"

  ## draw the layout of reads
  for i in xrange(min(MAX_READS_LAYOUT, len(reads))):
    r = reads[i]
    readid = r[0]
    start  = r[1]
    end    = r[2]
    rc     = r[3]
    color  = "blue"

    if (rc == 1): 
      color = "red"

    plt.plot ([start,end], [-2*i, -2*i], lw=4, color=color)

  ## draw the base of the coverage plot
  plt.plot([0, totallen], [0,0], color="black")

  ## draw the coverage profile
  for i in xrange(len(deltacov)-1):
    x1 = deltacov[i][0]
    x2 = deltacov[i+1][0]
    y1  = deltacov[i][1]
    y2  = deltacov[i+1][1]

    ## draw the horizonal line
    plt.plot([x1, x2], [y1, y1], color="black")
    
    ## and now the right vertical to the new coverage level
    plt.plot([x2, x2], [y1, y2], color="black")

  ## draw a vertical bar with the max coverage
  plt.plot([maxcovpos, maxcovpos], [2*maxcov, -(2*len(reads)+20)], color="green")

  plt.draw()
  plt.show()



## Summary
###############################################################################
##
## In this example, the heap-based plane sweep algorithm was 300 to 400 times
## faster than the basic brute force approach on my laptop! The advantage
## of the heap-based approach over the list based approach will become even 
## more significant with deeper coverage as log(D) and D become more significant
## differences.
##
## In the next session we will explore the related question of how to index 
## the alignments to be able to quickly determine which reads span a particular 
## location
##
###############################################################################
