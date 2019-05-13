For those curious, deepTools needs a new interval tree backend that support metadata associated with each interval. I previously made such a thing, called libGTF. Consequently, I'm just working on a (A) a python front-end for that and (B) some modifications specific to deepTools (namely, every interval needs an associated `deepTools_group` tag and exon bounds will be a new attribute associated with transcripts).

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.

ktring.h and kseq.h are from [Heng Li](http://lh3lh3.users.sourceforge.net/) and are available under an MIT license.

Usage
=====

The only class contained here is `GTF` and it only has only one function that should ever be used, `findOverlaps`.

Note that as is the case in deepTools, this package attempts to convert between chromosome naming systems. Because the conversion may not always be obvious, this can fail.

The GTF class
-------------

To read one or more files into an interval tree, one initializes a new `GTF` class:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF("some_file.gtf")

Multiple files can also be used:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF(["some_file.gtf", "some_other_file.bed.gz"])

Files may be optionally compressed and the compression magic number is used to determine this.

For GTF and BED12 files, exons are not stored by default, this can be changed with the `keepExons` option:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF(["some_file.gtf", "some_other_file.bed.gz"], keepExons=True)

The utility of this will be seen later. GTF and BED files may contain comments or browser lines at the beginning, these are ignored.

### Labels

It's often useful to have multiple groups of intervals. This can be accomplished by assigning a label to each interval. If multiple files are used, then this package will default to assigning the file name as a label to intervals in each input file. Alternatively, labels can be included inside of files. For BED files, this is accomplished as follows:

    chr1	1	100
    chr1	150	200
    #My group
    chr1	300	400
    #My other group

These labels **MUST** be unique in BED files. If they are not, then each subsequent instance will have a suffix appended to ensure that it is unique.

For GTF files, labels are included in the attribute column, by the addition of `deepTools_group` key:value pair:

    chr1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; deepTools_group "group 1";

These labels do **NOT** need to be unique across files.

Labels can be over-riden with the `labels` option:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF(["some_file.gtf", "some_other_file.bed.gz"], keepExons=True, labels=["foo", "bar", "quux", "sniggly"])

The number of provided labels **MUST** match the number encountered. These labels are applied in the order that groups are encountered in the input files. So if in the above example both files contain two groups, then the following would produce the same results but with the labels swapped across files:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF(["some_other_file.bed.gz", "some_file.gtf"], keepExons=True, labels=["foo", "bar", "quux", "sniggly"])

Labels can also be replaced after the fact:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF(["some_file.gtf", "some_other_file.bed.gz", "some_file.gtf"], keepExons=True)
    >>> gtf.labels = ["foo", "bar", "quux", "sniggly"]

### GTF-specific options

GTF files come with three options specific to them: `exonID`, `transcriptID`, and `transcript_id_designator`. The "feature" column (column 3) in a GTF file denotes the type of feature an entry describes. By default, this package only looks at entries with `transcript` or `exon` (with `keepExons=True`) in the feature column. For some use cases, one might instead want to store CDS as exonic intervals or replace transcripts with genes.

Transcripts are the primary entry used by this package and, consequently, each needs to have an associated transcript ID. Duplicate IDs are always ignored, since such a thing would be biologically non-sensical. In GTF files, the transcript id is stored in as `transcript_id "some ID";`. If, however, one changes thr `transcriptID` value to something else, such as `gene`, this key:value pair may not longer be present or may not be unique. In such cases, it's beneficial to change the key portion, for example to `gene_id`.

### Finding overlaps

Finding overlaps requires a chromosome, start, and end positions. As with BED files, these coordinates are 0-based half-open. By default, strand and overlap type are completely ignored. This can be overridden:

     >>> o = gtf.findOverlaps("chr1", 0, 100, strand="+", matchType=1, strandType=3)

This would search for intervals on the `+` strand (ignoring those on `.`, which would have additionally been returned had `strandType=1` been used) that are exactly [0, 100) on chromosome 1. Anyone interested in these more advanced overlap searching methods should look at the gtf.h file and the "libGTF" repository for examples.

It's often the case the a function looking for intervals is doing so by first dividing the genome into chunks and then sending each chunk to a processor for subsequent analysis. In such cases, it's convenient to NOT have processor duplicate processing intervals that may overlap multiple genomic bins. In these circumstances, the `trimOverlap` option can be set to `True`.

The output of `findOverlaps()` is a list of tuples:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF("foo.gtf", keepExons=True)
    gtf.findOverlaps("chr1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 12227), (12612, 12721), (13220, 14409)], '.'), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 12057), (12178, 12227), (12612, 12697), (12974, 13052), (13220, 13374), (13452, 13670)], '.'), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 14501), (15004, 15038), (15795, 15947), (16606, 16765), (16857, 17055), (17232, 17368), (17605, 17742), (17914, 18061), (18267, 18366), (24737, 24891), (29533, 29570)], '.'), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)], '.')]

Each tuple contains the following members (in order): 0-based starting position, 1-based end position, ID (the transcript ID for GTF files, column 4 for BED6/12 files and a string composed of the intervals for BED3 files), a group label, a sorted list of exonic bounds, and the score (column 4 in GTF files and 5 in BED files). If either the input file type does not provide exonic bounds or `keepExons=True` was not used, these bounds will be identical to that in the tuple:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF("foo.gtf")
    >>> gtf.findOverlaps("chr1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 14409)], '.'), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 13670)], '.'), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 29570)], '.'), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)], '.')]

In some cases, it's desirable to have the group labels be numeric, since the regions may be used for further processing and the results sorted or grouped accordingly. The `numericGroups` argument can be used to facilitate this:

    >>> from deeptoolsintervals import GTF
    >>> gtf = GTF("foo.gtf")
    >>> gtf.findOverlaps("chr1", 1, 20000, numericGroups=True)
    [(11868, 14409, 'ENST00000456328', 0, [(11868, 14409)], '.'), (12009, 13670, 'ENST00000450305', 0, [(12009, 13670)], '.'), (14403, 29570, 'ENST00000488147', 0, [(14403, 29570)], '.'), (17368, 17436, 'ENST00000619216', 1, [(17368, 17436)], '.')]

The Enrichment class
--------------------

The `Enrichment` class is a modification of the base `GTF` class, aimed at querying feature types in a region. Creation of the class from one or more BED/GTF files is also similar to the `GTF` class:

    >>> from deeptoolsintervals import Enrichment
    >>> gtf = Enrichment(["foo.gtf", "bar.bed"])

For GTF files, the feature type is the 3rd column. For BED files, the feature type is the file name, though this can be changed with the `labels` option:

    >>> from deeptoolsintervals import Enrichment
    >>> gtf = Enrichment(["foo.gtf", "bar.bed"], labels=["this will be ignored", "peaks"])

For GTF files, the label is ignored, but for the sake of simplicity if labels are specified there must be at least as many as there are files. Note that, unlike with `GTF` objects, you can not change labels after creation of an `Enrichment` object.

All entries in BED and GTF files are stored. For BED12 files, only columns 2/3 are used as region bounds by default. This can be modified with the `keepExons` option:

    >>> from deeptoolsintervals import Enrichment
    >>> gtf = Enrichment("bar.bed12", keepExons=True)

All other file types will ignore the `keepExons` option, so this can still be specified with a mix of BED12 and other file types.

### Finding overlaps

Finding overlaps of an `Enrichment` object is similar to that with a `GTF` object. Once again, the `findOverlaps()` method is used, though the `trimOverlap`, `numericGroup`, and `includeStrand` options are not present. Further, instead of a single `start` and `end` value, a list of tuples is used. This last difference facilitates finding overlaps of spliced genes, since the pysam `get_blocks()` method returns this type of data:

    >>> from deeptoolsintervals import Enrichment
    >>> gtf = Enrichment("GRCh38.84.gtf.gz")
    >>> gtf.findOverlaps("1", [(65500, 65600), (69900, 70000)])
    frozenset(['start_codon', 'transcript', 'gene', 'exon', 'CDS'])

The output is a set containing all overlapped feature types. This is convenient for quick summarization.

## Enrichment of custom attributes

As of deeptoolsintervals 0.1.8, the `Enrichment` class is able to use a custom attribute key instead of the feature type. This allows you to find overlaps of things like the gene biotype:

    >>> from deeptoolsintervals import Enrichment
    >>> gtf = Enrichment("GRCh38.84.gtf.gz", keepExons=True, attributeKey="gene_biotype")
    >>> gtf.findOverlaps("1", [(0, 2000000)])
    frozenset(['miRNA', 'group 1', 'group 2', 'transcribed_unprocessed_pseudogene', 'processed_pseudogene', 'lincRNA', 'unprocessed_pseudogene', 'protein_coding']))
