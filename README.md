# DC-TAP-seq Internal Toolkit

> [!NOTE]
> This repo is currently undergoing development.
> Most of the scripts and analyses currently committed are outdated. As such,
> a restructuring and refactoring effort of this codebase is under way.

Hi! Thanks for coming to visit this repo!
My intent with this repository is to curate a collection of tools that enables
the design of large-scale, CRISPRi-based, perturbation screening experiments
similar to one reported in the DC-TAP-seq paper
(Ray J., _et al. Biorxiv_ 2025.09.16.676677.
DOI: [10.1101/2025.09.16.676677](https://doi.org/10.1101/2025.09.16.676677)).

This project focuses on organizing existing tools and creating new ones that are
both flexible and robust. Each tool should be usable on its own or as
part of a uniform pipeline or workflow. Users should be able to design large-
scale CRISPRi guide libraries, optimize primer design, run basic validation
analyses, and generate standard QC plots that enable quick, intuitive assessment
of data quality for decision making. These plots maybe helpful in indicating
whether generated datasets are strong enough for expert analysis,
need further work (and possible rescue), or are too poor to continue.
My hope for these tools is to get any interested users at least 80% of the
way towards designing, validating, and analyzing large-scale perturbation
experimental datasets.

## Acknowledgements and Disclaimers

Please note that all credit for referenced works belongs to the original
authors and their institutions. See the official
[EngreitzLab/DC_TAP_Paper repo](https://github.com/EngreitzLab/DC_TAP_Paper),
published papers, and
[preprints](https://www.biorxiv.org/content/10.1101/2025.09.16.676677v1)
for upstream work that has motivated the efforts here.

**This work is largely a personal effort.** Most of the code and analyses have
not been formally reviewed or validated by collaborators, organizations, or
institutions. While I plan to rigorously refine these tools for future use
(personal, career, and institutional), please be mindful of the non-standard
validations and _use_ them _at your own discretion_.

That said—I want to take a moment to thank the many mentors, collaborators,
supporters, and friends who have guided me on this journey. A huge shout out
to the folks at the
[Broad Institute of MIT and Harvard](https://www.broadinstitute.org/),
the HGRM Enhancer Team at the
[Novo Nordisk Foundation Center for Genomic Mechanisms of Disease](https://www.broadinstitute.org/nnfc),
[Jesse Engreitz](https://med.stanford.edu/profiles/jesse-engreitz) and his
[lab at Stanford](https://www.engreitzlab.org/), the Shoresh computational group
under the Broad’s Epigenomics Program, and the Jones and Donnard groups at the
[Lander Lab](https://www.broadinstitute.org/lander-lab/people).
And last but definitely not least: a heartfelt thanks to
[Judhajeet Ray](https://www.broadinstitute.org/bios/judhajeet-ray), whose
guidance has pushed me to always be the best scientist I can and should be.

None of this work would have been possible without such generous mentorship
and support.
