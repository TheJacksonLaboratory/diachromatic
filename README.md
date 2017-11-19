# diachromatic
Differential analysis of chromatin interactions by capture

## Summary
This package will implement parts of the capture Hi-C pipeline as done in
HiCUP, but will do so in a way that is easily integrated with VPV. The final
part of the pipeline will perform differential analysis of 
promoter reads using a likelihood ratio test.

### Digest
The first part of the pipeline creates a digest file representing
the in silico digestion of the genome using the same enzyme(s) that
the user chose in VPV for the capture panel.

```
$ java -jar Diachromatic.jar digest -g genomeDir [-o filename]
```
When run in the digest mode, Diachromatic needs to be passed the
path to the directory that contains the genome FASTA files. This
should be the same directory used to choose the viewpoints with VPV.
The default name for the output file is XXX, but can be overridden with 
the -o option.
