Pretrained Q matrices:
- Q.mammal: trained on mammalian dataset realigned with muscle v5.1 ([Edgar 2022])
- Q.mammal_cloak: trained on mammalian dataset realigned with muscle v5.1, and filtered using cloak
- Q.mammal_guidance: trained on mammalian dataset realigned with muscle v5.1, and filtered using GUIDANCE v2.01 ([Sela, et al. 2015]) by masking all residues with residue scores < 0.9 
- Q.mammal_taper: trained on mammalian dataset realigned with muscle v5.1, and filtered using Taper ([Zhang, et al. 2021])
- Q.mammal_divvier: trained on mammalian dataset realigned with muscle v5.1, and filtered using divvier ([Ali, et al. 2019]) with the -divvygap option
- QC.mammal: trained on mammalian dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap and -partial options
- Q.mammal_GTRPMIX: trained on mammalian dataset realigned with muscle v5.1 under a +C60 profile mixture model with linked exchangeabilities ([Banos, et al. 2024]). It is optimized to be used for profile mixture models
- QC.arch: trained on archaeal dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap and -partial options
- QC.bac: trained on bacterial dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap and -partial options
- QC.pfam: trained on pfam dataset realigned with muscle v5.1, and filtered using divvier with the -divvygap and -partial options

The training datasets for these matrices can all be found in the [Data](/Data) folder of this repository
All QC matrices are also available on IQTree

[Edgar 2022]: https://doi.org/10.1038/s41467-022-34630-w
[Sela, et al. 2015]: https://doi.org/10.1093/nar/gkv318
[Zhang, et al. 2021]: https://doi.org/10.1111/2041-210X.13696
[Ali, et al. 2019]: https://doi.org/10.1093/molbev/msz142
[Banos, et al. 2024]: https://doi.org/10.1093/molbev/msae174
