h1. Genomic Peptide Finder

The Genomic Peptide Finder is a tool that matches de novo predicted amino 
acid sequences to the genomic DNA sequence of an organisms. The matching
procedure is error-tolerant to accomodate possible sequencing errors due
to individual missing fragment ions. In addition, spliced peptides can be
deduced.

GPF has been devised by Jens Allmer and Michael Hippler in 2003, with 
subsequent design changes and reimplementation by Michael Specht in 2007.


h2. Building GPF

GPF requires Qt 4.6+, and building is straightforward:

<pre>
$ cd projects/gpfindex
$ qmake
$ make release

$ cd projects/gpfquery
$ qmake
$ make release
</pre>

