#!/usr/bin/env python

from deeptoolsintervals import tree
import re
import sys
import gzip

gene_id_regex = re.compile('(?:gene_id (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])')
transcript_id_regex = re.compile('(?:transcript_id (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])')
deepTools_group_regex = re.compile('(?:deepTools_group (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])')


def seemsLikeGTF(cols):
    """
    Does a line look like it could be from a GTF file? Column contents must be:

    3: int
    4: int
    5: '.' or float
    6: '+', '-', or '.'
    7: 0, 1 or 2
    8: matches the attribute regular expression
    """
    try:
        int(cols[3])
        int(cols[4])
        if cols[5] != '.':
            float(cols[5])
        cols[6] in ['+', '-', '.']
        if cols[7] != '.':
            int(cols[7]) in [0, 1, 2]
        assert(gene_id_regex.match(cols[8]) != None)
        return True
    except:
        return False


def findRandomLabel(labels, name):
    """
    Because some people are too clever by half, ensure that group labels are unique...
    """
    if name not in labels:
        return name

    # This is what the heatmapper.py did to ensure unique names
    i = 0
    while True:
        i += 1
        name = name + "_r" + str(i)
        if name not in labels:
            return name


def parseExonBounds(start, end, n, sizes, offsets):
    """
    Parse the last 2 columns of a BED12 file and return a list of tuples with
    (exon start, exon end) entries.

    If the line is malformed, issue a warning and return (start, end)
    """
    offsets = offsets.strip(",").split(",")
    sizes = sizes.strip(",").split(",")
    try:
        starts = [start + int(x) for x in offsets]
        ends = [start + int(x) + int(y) for x, y in zip(offsets, sizes)]
    except:
        sys.stderr.write("Warning: Received an invalid exon offset ({}) or size ({}), using the entry bounds instead ({}-{})\n".format(offsets, sizes, start, end))
        return [(start, end)]

    if len(offsets) < n or len(sizes) < n:
        sys.stderr.write("Warning: There were too few exon start/end offsets ({}) or sizes ({}), using the entry bounds instead ({}-{})\n".format(offsets, sizes, start, end))
        return [(start, end)]

    return [(x, y) for x, y in zip(starts, ends)]


class GTF(object):
    """
    A class to hold an interval tree and its associated functions
    """

    def __init__(self):
        """
        fnames:	The file name(s)
        exonID:	The exon feature designator (only for GTF)
        transcriptID: The transcript feature designator (only for GTF)
        labels:	A list of group labels (designated by 'deepTools_group' in GTF or # lines in BED or multiple file names for both)
        ftype:	The (current) inferred file type
        chroms:	A dictionary of chromosomes that have been found. This is only needed in case someone uses multiple files with different chromosome names
        exons:	A dictionary with transcript ID as the key and exonic bounds as the value. We'll see if this is quick enough
        labelList: The input labels, if they exist
        tree: The pyGTFtree object is held here
        """
        self.fnames = []
        self.exonID = "exon"
        self.transcriptID = "transcript"
        self.keepExons = False
        self.chroms = []
        self.exons = {}
        self.labels = []
        self.labelList = []
        self.transcriptIDduplicated = []
        self.tree = tree.initTree()
        self.labelIdx = 0

    def firstNonComment(self, fp):
        line = fp.next()
        try:
            while line.startswith("#") or line.startswith('track') or line.startswith('browser'):
                line = fp.next()
        except:
            sys.stderr.write("Warning, {} was empty\n".format(fp.name))
            return None
        return line

    def inferType(self, fp, line):
        cols = line.split("\t")
        if len(cols) < 3:
            raise RuntimeError('{} does not seem to be a recognized file type!'.format(fp.name))
        elif len(cols) == 3:
            return 'BED3'
        elif len(cols) < 6:
            sys.stderr.write("Warning, {} has an abnormal number of fields. Assuming BED3 format.\n".format(fp.name))
            return 'BED3'
        elif len(cols) == 6:
            return 'BED3'
        elif len(cols) == 9 and seemsLikeGTF(cols):
            return 'GTF'
        elif len(cols) == 12:
            return 'BED12'
        elif len(cols) < 12:
            sys.stderr.write("Warning, {} has an abnormal format. Assuming BED6 format.\n".format(fp.name))
            return 'BED6'
        else:
            sys.stderr.write("Warning, {} has an abnormal format. Assuming BED12 format.\n".format(fp.name))
            return 'BED12'

    def mungeChromosome(self, chrom):
        """
        Return the chromosome name, possibly munged to match one already found in the chromosome dictionary
        """
        if chrom in self.chroms:
            return chrom

        # chrM <-> MT and chr1 <-> 1 conversions
        if chrom == "MT" and "chrM" in self.chroms:
            chrom = "chrM"
        elif chrom == "chrM" and "MT" in self.chroms:
            chrom = "MT"
        elif chrom.startswith("chr") and len(chrom) > 3 and chrom[3:] in self.chroms:
            chrom = chrom[3:]
        elif "chr" + chrom in self.chroms:
            chrom = "chr" + chrom

        self.chroms.append(chrom)

        return chrom

    def parseBEDCore(self, line, ncols):
        """
        Returns True if the entry was added, otherwise False
        """

        strand = 3
        cols = line.split("\t")
        name = "{}:{}-{}".format(cols[1:3])

        if int(cols[1]) >= int(cols[2]) or int(cols[1]) < 0:
            sys.stderr.write("Warning: {}:{}-{} is an invalid BED interval! Ignoring it.\n".format(cols[0], cols[1], cols[2]))
            return

        # BED6/BED12: set name and strand
        if ncols != 3: 
            name = cols[4]
            if cols[6] == '+':
                strand = 0
            elif cols[6] == '-':
                strand = 1

        # Ensure that the name is unique
        if name in self.exons.keys():
            sys.stderr.write("Skipping {}, an entry by this name already exists!\n".format(name))
            return False
        else:
            self.tree.addEntry(self.mungeChromosome(cols[0]), int(cols[1]), int(cols[2]), name, strand, self.labelIdx)
            if ncols != 12 or self.keepExons == False:
                self.exons[name] = [(int(cols[1]), int(cols[2]))]
            else:
                self.exons[name] = parseExonBounds(int(cols[1]), int(cols[2]), int(cols[9]), cols[10], cols[11])
        return True

    def parseBED(self, fp, line, ncols=3):
        """
        parse a BED file. The default group label is the file name.
        """
        groupLabelsFound = 0
        groupEntries = 0
        # If are already labels then we need to increment the index, otherwise we start overwritting
        # TODO What happens if we use multiple files, one of which is only comments?
        if self.labelIdx > 0:
            self.labelIdx += 1
        startingIdx = self.labelIdx

        # Handle the first line
        if parseBEDcore(line, ncols):
            groupEntries = 1

        # iterate over the remaining lines
        for line in fp:
            if line.startswith("#"):
                # If there was a previous group AND it had no entries then remove it
                if groupLabelsFound > 0:
                    if groupEntries == 0:
                        sys.stderr.write("Warning, the '{}' group had no valid entries! Removing it.\n".format(self.labels[self.labelIdx]))
                        del self.labels[-1]
                        groupLabelsFound -= 1
                        self.labelIdx -= 1

                label = line[1:].strip()
                if len(label):
                    # Gaard against duplicate group labels
                    self.labels.append(findRandomLabel(self.labels, label))
                else:
                    # I'm sure someone will try an empty label...
                    self.labels.append(findRandomLabel(self.labels, fp.name))
                self.labelIdx += 1
                groupLabelsFound += 1
                groupEntries = 0
            else:
                if parseBEDcore(line, ncols):
                    groupEntries += 1

        if groupLabelsFound == 0 or self.labelIdx - startingIdx > groupLabelsFound:
            # This can only happen once
            self.labels.append(findRandomLabel(self.labels, fp.name))
            self.labelIdx += 1

    def parseGTFtranscript(self, cols, label):
        """
        Parse and add a transcript entry
        """
        if int(cols[3]) - 1 < 0:
            sys.stderr.write("Warning: Invalid start in '{}', skipping\n".format("\t".join(cols)))
            return

        if len(cols) < 9:
            sys.stderr.write("Warning: non-GTF line encountered! {}\n".format("\t".join(cols)))
            return

        m = deepTools_group_regex.search(cols[8])
        if m:
            label = m.group(1)

        m = transcript_id_regex.search(cols[8])
        if not m:
             sys.stderr.write("Warning: {} is malformed!\n".format("\t".join(cols)))
             return

        name = m.group(1)
        if name in self.exons.keys():
            sys.stderr.write("Warning: {} occurs more than once! Only using the first instance.\n".format(name))
            self.transcriptIDduplicated.append(name)
            return

        if int(cols[3]) > int(cols[4]) or int(cols[3]) < 1:
            sys.stderr.write("Warning: {}:{}-{} is an invalid GTF interval! Ignoring it.\n".format(cols[0], cols[3], cols[4]))
            return

        strand = 3
        if cols[6] == '+':
            strand = 0
        elif cols[6] == '-':
            strand = 1

        # Get the label index
        if label not in self.labels:
            self.labels.append(label)
        self.labelIdx = self.labels.index(label)

        chrom = self.mungeChromosome(cols[0])
        self.tree.addEntry(chrom, long(cols[3]) - 1, long(cols[4]), name, long(strand), long(self.labelIdx))

        # Exon bounds placeholder
        self.exons[name] = []

    def parseGTFexon(self, cols):
        """
        Parse an exon entry and add it to the transcript hash
        """
        if int(cols[3]) - 1 < 0:
            sys.stderr.write("Warning: Invalid start in '{}', skipping\n".format("\t".join(cols)))
            return

        m = transcriptID_regex.search(cols[9])
        if not m:
            sys.stderr.write("Warning: {} is malformed!\n".format("\t".join(cols)))
            return

        name = m.group(0)
        if name in self.transcriptIDduplicated.keys():
            return
        if name not in self.exons.keys():
            self.exons[name] = []

        self.exons[name].append((int(cols[3]) - 1, int(cols[4])))

    def parseGTF(self, fp, line):
        """
        parse a GTF file. Note that a single label will be used for every entry
        in a file that isn't explicitly labeled with a deepTools_group
        key:values pair in the last column
        """
        file_label = findRandomLabel(self.labels, fp.name)

        # Handle the first line
        cols = line.split("\t")
        if cols[2].lower() == self.transcriptID:
            self.parseGTFtranscript(cols, file_label)
        elif cols[2].lower() == self.exonID:
            self.parseGTFexon(cols)

        # Handle the remaining lines
        for line in fp:
            if not line.startswith('#'):
                cols = line.split("\t")
                if cols[2].lower() == self.transcriptID:
                    self.parseGTFtranscript(cols, file_label)
                elif cols[2].lower() == self.exonID and self.keepExons == True:
                    self.parseGTFexon(cols)

        # Reset self.labelIdx
        self.labelIdx = len(self.labels) - 1

    def parse(self, fnames, exonID="exon", transcriptID="transcript", keepExons=False, labelList=[]):
        """
        Driver function to actually parse files. The steps are as follows:

        1) skip to the first non-comment line
        2) Infer the type from that
        3) Call a type-specific processing function accordingly
          * These call the underlying C code for storage
          * These handle chromsome name conversions (python-level)
          * These handle labels (python-level, with a C-level numeric attribute)
        4) Sanity checking (do the number of labels make sense?)
        """
        self.exonID = exonID
        self.transcriptID = transcriptID
        self.keepExons = keepExons
        self.labels = labelList

        if labelList != []:
            self.already_input_labels = True

        # Load the files
        for fname in fnames:
            fp = gzip.open(fname, "rb")
            line = self.firstNonComment(fp)
            if line is None:
                # This will only ever happen if a file is empty or just has a header/comment
                continue

            ftype = self.inferType(fp, line)
            if ftype == 'GTF':
                self.parseGTF(fp, line)
            elif ftype == 'BED':
                self.parseBED(fp, line, 3)
            elif ftype == 'BED6':
                self.parseBED(fp, line, 6)
            else:
                self.parseBED(fp, line, 12)
            fp.close()

        # TODO Sanity check
        # TODO Are there entries?
        # TODO Do the number of groups match the number of labels?
        # Replace labels if there's a labelList

        # vine -> tree
        self.tree.finish()

   # findOverlaps()
