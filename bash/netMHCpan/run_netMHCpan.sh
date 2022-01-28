
## netMHCpan 4.0
cd /Users/janihuuh/Dropbox/applications/netMHCpan-4.0/run/

## Download the FASTA from uniprot.org, e.g. https://www.uniprot.org/uniprot/Q16655.fasta
../netMHCpan mlana.fsa -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01 > mlana.fsa.myout;
../netMHCpan pp65.fsa -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01 > pp65.fsa.myout;
../netMHCpan meloe1.fsa -a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01 > meloe1.fsa.myout;
