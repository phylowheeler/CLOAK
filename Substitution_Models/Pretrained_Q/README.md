Pretrained Q matrices:
- Q.mammal: trained on mammalian dataset realigned with muscle v5.1
- Q.mammal_cloak: trained on mammalian dataset realigned with muscle v5.1, and filtered using cloak
- Q.mammal_guidance: trained on mammalian dataset realigned with muscle v5.1, and filtered using GUIDANCE v2.01 by masking all residues with residue scores < 0.9 ([Sela, et al. 2015])
- Q.mammal_taper: trained on mammalian dataset realigned with muscle v5.1, and filtered using Taper
- Q.mammal_divvier: trained on mammalian dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap option
- QC.mammal: trained on mammalian dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap and -partial options
- Q.mammal_GTRPMIX: trained on mammalian dataset realigned with muscle v5.1 under a +C60 profile mixture model with linked exchangeabilities. It is optimized to be used for profile mixture models
- QC.arch: trained on archaeal dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap option
- QC.bac: trained on bacterial dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap option
- QC.pfam: trained on pfam dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap option

All QC matrices are also available on IQTree

[Sela, et al. 2015]: https://doi.org/10.1093/nar/gkv318
