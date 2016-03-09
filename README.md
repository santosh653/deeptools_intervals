For the love of <insert diety> don't use this repository! It only exists to more conveniently create a python interface to libGTF for deepTools!
================================================================================================================================================

For those curious, deepTools needs a new interval tree backend that support metadata associated with each interval. I previously made such a thing, called libGTF. Consequently, I'm just working on a (A) a python front-end for that and (B) some modifications specific to deepTools (namely, every interval needs an associated `deepTools_group` tag and exon bounds will be a new attribute associated with transcripts).

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.

ktring.h and kseq.h are from [Heng Li](http://lh3lh3.users.sourceforge.net/) and are available under an MIT license.
