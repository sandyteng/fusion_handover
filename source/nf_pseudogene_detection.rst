=======================================
Pseudogene detection pipeline handover
=======================================

-----------------
Purpose
-----------------

This document records pipeline description, installation and execution

----

-----------------
Design Concept
-----------------

This workflow/pipeline employs three aligners—BLAT, MegaBLAST, and BLASTN—to perform a two-stage pseudogene detection process. The input consists of amplicon sequences suspected to contain pseudogene regions (due to non-specific alignments). Any such sequences identified will be recorded in the output file pseudolist.report, which lists all amplicons affected by potential pseudogenes.

In the first stage, the pipeline detects amplicons that produce single-hit alignments using both MegaBLAST and BLAT. (A single-hit amplicon is defined as having only one alignment result from both aligners.) Amplicons with multiple alignments in either MegaBLAST or BLAT proceed to the second stage.

In the second stage, these ambiguous cases are further analyzed using BLASTN. If BLASTN reports two or more alignments for a given read or amplicon, it is classified as pseudogene-affected and included in pseudolist.report.

----

-----------------
Repo
-----------------

- `Pseudogene_detection <https://github.com/ACTGenomics/Pseudogene_detection>`_


-----------------
Docker Image
-----------------

Stable version for Pseudogene detection pipeline: `sandyteng/pseudogene_sandy:v0.0.1 <https://hub.docker.com/repository/docker/sandyteng/pseudogene_sandy/general>`_
Backup image (identical to 'sandyteng/pseudogene_sandy:v0.0.1'): `actgenomics/pseudogene_sandy:v0.0.1 <https://hub.docker.com/repository/docker/actgenomics/pseudogene_sandy/general>`

Remark:
If the 'sandyteng/pseudogene_sandy:v0.0.1' is not available, one can use the 'actgenomics/pseudogene_sandy:v0.0.1' instead (Modify 'container' in the configuration file).

.. code-block:: console

    container = 'docker.io/sandyteng/pseudogene_sandy:v1.0.0'

----

-----------------
Others materials
-----------------

- Design document - Pseudogene detection: `pseudogene identification.pptx <https://actgenomics-my.sharepoint.com/:p:/p/sandyteng/EbqdP70a7EVHlz_HK0YIStIBxZMUdUx3thEsO87q4qTk_w?e=D8UMBH>`_

----

--------------------
Execution
--------------------

Server should already have native nextflow installed, if not uses conda to create a Nextflow run environment

.. code-block:: console

    ### The following commands are copied from the latest testing runs (ref. issue: https://actg.atlassian.net/browse/ABIE-836)

    # Execute the pipeline to classify PA039(GRCh38) amplicon sequences
    nextflow run /mnt/RD_Develop/sandyteng/workdir/repo_test/Pseudogene_detection/main_grch38.nf \
        -params-file /mnt/RD_Develop/sandyteng/workdir/repo_test/testparams/PA039_GRCh38.json \
        -c /mnt/RD_Develop/sandyteng/workdir/repo_test/Pseudogene_detection/pseudogene_localdocker.config

    # To-do (PA031(hg19), Onco2M7(hg19))

--------------------
Conclusion
--------------------

This completes the instructions for running Fusion V5 pipelines.