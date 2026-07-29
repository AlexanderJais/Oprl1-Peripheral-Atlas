# Figure legends

## Figure 1. *Oprl1* is the highest-expressed opioid receptor gene in peripheral neurons

Mean expression of the four opioid receptor genes in each of 19 peripheral neuronal populations,
grouped by division of the peripheral nervous system. Bars are pseudobulk means over all neurons
of the population, in the unit each dataset was measured in, ordered within each panel from
highest to lowest. *Oprl1* is blue and the other three receptors are gray. Sample size is the
number of neurons contributing to the panel.

(A) Spiral ganglion, GSE114997, full-length SMART-seq, median 2,678,701 reads per neuron. *Oprl1*
41.94 CPM in 74.8% of neurons; *Oprm1*, *Oprd1*, and *Oprk1* detected in no cell. All four genes
are present in the annotation, so the three zeros are measured absences.

(B and D) Geniculate ganglion on two platforms from two laboratories. (B) GSE102443, full-length
SMART-seq, FPKM, 31.30-fold over the second receptor, bootstrap support 1.000. (D) GSE135801, 3'
droplet, 4.02-fold, support 0.974 (95% interval 0.99 to 36.70).

(C) Vestibular ganglion, GSE309608, four mice, 4.29-fold (3.85 to 4.81), support 1.000, same
ordering in each mouse.

(E) Trigeminal ganglion, iPain atlas, whole-cell neurons only, 2.08-fold (1.41 to 3.26).

(F and G) Nodose and jugular ganglia, NodoMap atlas, whole-cell datasets only. Nodose 1.21-fold
(1.15 to 1.28), support 1.000, first in all four contributing datasets; jugular 1.16-fold (0.92 to
1.46), support 0.898, first in two of three.

(H) Dorsal root ganglion, iPain atlas, whole-cell neurons only, 1.13-fold (1.08 to 1.19), the
narrowest margin measured and the second largest sample.

(I to L) Sympathetic ganglia. Celiac 44.41-fold (20.42 to 197.63), superior cervical 25.72-fold
(15.93 to 52.14), thoracic chain 14.69-fold (9.35 to 27.27), stellate 4.68-fold (3.02 to 8.08).
Support 1.000 throughout. (I) is from GSE232789, (J) from the two untreated animals of GSE231766,
(K) from GSE78845, and (L) from GSE231924, which reports cardiac-projecting neurons and places
*Oprl1* first in each of eight mice.

(M and N) Parasympathetic ganglia. (M) Sphenopalatine ganglion, GSE232789, 66.86-fold (40.54 to
132.54), the largest margin in the study. These neurons are cholinergic and not noradrenergic
(*Slc18a3* 252 CPM, *Chat* 26, *Th* 2, *Dbh* 42) against *Th* 390 to 840 CPM in the sympathetic
ganglia of the same experiment. (N) Intrinsic cardiac nervous system, GSE330884, three mice,
3.65-fold (3.29 to 4.08), the only population in this division in which *Oprm1* is the second
receptor.

(O and P) Enteric submucosal neurons of the small intestine at two ages, GSE263422, one laboratory
and one platform. (O) Postnatal day 24, *Oprl1* first at 2.20-fold (2.05 to 2.37) in all three
samples. (P) Postnatal day 7, *Oprk1* first at 122.51 CPM against *Oprl1* 26.70, 4.59-fold (4.04
to 5.24) in both samples. *Oprk1* falls to 3.77 CPM at 1.1% detection by day 24 while *Oprl1*
holds between 26 and 36 CPM.

All panels use whole-cell data. Both iPain atlases are majority single-nucleus and are restricted
here to `suspension_type == "cell"`; pooling preparations compresses the trigeminal margin from
2.08-fold to 1.18-fold, for the reason given in Figure S1. Neurons were called on raw counts with
glial and immune barcodes excluded, and every population passed the same marker gate and
ambient-RNA check before any receptor value was read from it. Margins are the ratio of the highest
to the second-highest receptor; intervals in parentheses are 95% bootstrap intervals over 10,000
resamples of the cells, and support is the fraction of those resamples retaining the observed top
receptor. Five populations (I, M, and the lumbar chain, pelvic, and stellate ganglia of Table S1)
come from a single experiment in which six autonomic ganglia were dissected and sequenced
together.

See also Figure S1, Figure S3, and Table S1.
