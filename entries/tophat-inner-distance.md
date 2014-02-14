A common way of estimating the ``--mate-inner-dist`` option for Tophat is to align a small subset of reads with Bowtie and calculate the insert size based on the TLEN field. This is a bit of a fudge because some of the read pairs cross exon-exon boundaries and those reads will have a very large value for TLEN. One way to get around this is to filter out unreasonably large values of TLEN but that assumes you know ahead of time how large the insert sizes should be.

Using the median instead of the mean and an estimate of the standard deviation based on the median absolute deviation solves is more robust to the exon-exon outliers. Here is some example code to perform this calculation:


```
from contextlib import closing
import pysam
import numpy

def innerdist(sam_file):
    dists = []
    with closing(pysam.Samfile(out_sam)) as work_sam:
        for read in work_sam:
            if read.is_proper_pair and read.is_read1:
                dists.append(abs(read.isize) - 2 * read.len)
    median = float(numpy.median(dists))
    deviations = [abs(d - median) for d in dists]
    # median absolute deviation estimator of the standard deviation
    mad = 1.4826 * float(numpy.median(deviations))
    return int(median), int(mad)
```

Running ``innerdist`` on a final library with this Bioanalyzer plot:
![](http://dl.dropboxusercontent.com/u/2822886/Screenshots/04.png)  and 100 basepair reads returns values of -7 and 72 for the mean and standard deviation respectively, which pretty much nails it. The Illumina adapter sequences are about 120 bases in total, which means the insert being sequenced is 200 base pairs and with 100 basepair reads that means the inner distance should be around 0.

Using the mean and standard deviation returns estimates of these parameters as 79 and 255.
