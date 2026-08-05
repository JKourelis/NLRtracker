# NLRtracker 1.x changelog

Preserved from the original README. For 2.0.0 and later see [`../README.md`](../README.md).

## v1.0.3
- Update to InterProScan 5.53-87.0
- Improved FASTA read-in (using `Biostrings::readAAStringSet`)
- Improved identification of the C-JID domain (addition of the PF20160 signature)
- Extraction of sequences containing a C-JID domain but lacking an NB-ARC domain
- Corrected a mistake fusing multiple sequential NB-ARC domains
- Corrected assignment of CCG10-NLRs as CC-NLRs in the subclass column

## v1.0.2
- C-JID domain annotation using HMMER
- CblN domain (overlapping with RPW8-type CC domains)
- Output for use with iTOL

## v1.0.1
- Update to InterProScan 5.51-85.0
- Improved NB-ARC extraction (including the winged-helix domain)
- Improved identification of Rx-type CC (addition of the G3DSA:1.20.5.4130 signature)

## v1.0.0
- Initial NLRtracker release
