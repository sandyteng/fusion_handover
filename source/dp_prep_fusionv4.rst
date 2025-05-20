==============================================
Fusion V4-based workflow db files preparation
==============================================

-----------------
Purpose
-----------------

This document records how the Fusion V4 annotation files are prepared

----

-----------------
Repo
-----------------

- `Docker image construction for the commonly used fusion dev. tools <https://github.com/ACTGenomics/fusion_pipeline_env>`_

-----------------
Docker Image
-----------------

Stable version for the commonly used fusion dev. tools: `v0.6 <https://hub.docker.com/repository/docker/actgenomics/fusion_dev/general>`_

-----------------
Others materials
-----------------

- DB preparation summary (MANE v1.4): `[Fusion v5] preferred transcriptome fusion v4 -> v5 db update (MANE v0.95 -> v1.4) <https://actg.atlassian.net/browse/ABIE-1012>`_
- See the detailed db-v3.1 documentation here: `My PDF Guide <_static/Fusion_db_prep.steps_db-v3.1.pdf>`_
- See the executed commands for db-v3.1 transcriptome construction: `Executed commands for DB preparation steps <_static/Fusion_db_ref_transcript_v5_draft.pdf>`_

----

-----------------
Source Files
-----------------

These files are the source (input) files for fusion V4 db construction. 

- **MANE v1.4 DB**: 
    - ``MANE.GRCh38.v1.4.summary.txt.gz``
    - ``MANE.GRCh38.v1.4.ensembl_genomic.gff.gz``

- **Genome sequence, Grch38, GENCODE-r47**
    - ``GRCh38.p14.genome.fa.gz``

- **Probe information file provided by AD team**
    - ``ACTFusionv5_target-region_PartAB_individual_1039.bed``

- **Kinase files (manually curated)**
    - ``protein.26.v1.4.kinase.fasta``
    - ``protein.26.v1.4.kinase.meta.txt``

- **White list (gsp pair/probe pair)**
    - ``gsppairs_inclusion_v1.4.txt``
    - ``filter_internal.QC9.0.mgsp.qcr.0.5.blank.config``

.. code-block:: console

    # Local file list (directories)
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gff.gz

    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/gencode_v47/OpenDB_GENCODE_human_r47/GRCh38.p14.genome.fa.gz

    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/captureprobe_250401/ACTFusionv5_target-region_PartAB_individual_1039.bed

    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/v1.4_inputfiles/protein.26.v1.4.kinase.fasta
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/v1.4_inputfiles/protein.26.v1.4.kinase.meta.txt

    /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/data/gsppairs_inclusion_v1.4.txt
    /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/testconfigs/filter_internal.QC9.0.mgsp.qcr.0.5.blank.config

----

Workflows
~~~~~~~~~~~~~~
.. image:: _img/DB-v3.1_steps_1-4.png
    :width: 600px
    :align: center
    :alt: DB preparation part I

.. image:: _img/DB-v3.1_steps_5-8.png
    :width: 600px
    :align: center
    :alt: DB preparation part II

.. image:: _img/DB-v3.1_steps_9-11.png
    :width: 600px
    :align: center
    :alt: DB preparation part III
    
-----

-----------------
Tools
-----------------

- **Tool summary**
    - * indicates that some of the tools are available within image 'actgenomics/fusion_dev:v0.6'

.. code-block:: console

    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | No.    | Image | Tool(s)                                                                                                           |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 1      | *     | wget, rsync, zcat, samtools                                                                                       |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 2      | *     | awk, cat                                                                                                          |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 3      |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/filter_mane_gff.py                                                     |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 4      | *     | /tools/Fusion/convert2bed                                                                                         |
    |        |       | bedtools getfasta                                                                                                 |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 5(1)   |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/RefFusion.v2.py                                                        |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 8-0-a  |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/candidate_exons_mapping.sh                                            |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 8-0-b  |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/Probe_faheader_converter.py                                            |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 8-0-c  | *     | /tools/Fusion/ncbi-blast/bin/blastn                                                                               |
    |        |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/blastnparser.py                                                        |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 9      | *     | /tools/Fusion/bwa index                                                                                           |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+
    | 11     |       | /mnt/RD_Develop/sandyteng/ACTFusionV5/code/Get_shifted_boundary.py                                                |
    |        |       | /mnt/BI3/Team_workdir/sandyteng_workdir/ACTFusionV4_Torrent/code/update_qcconfig_with_tsv.py                     |
    +--------+-------+--------------------------------------------------------------------------------------------------------------------+

----

-----------------
Config Files
-----------------
- **Config file for v1.4 MANE Select transcriptome**
    - ``fusion_multi_localdocker.v9.20241125.v0.23.0_v1.4.MANE.transcriptome.v3-1.config``

- **Config file for v0.95 MANE Select transcriptome (== fusion v5 pipeline v0.1 config file for fusion v4-based workflow)**
    - ``fusion_multi_localdocker.v9.20241125.v0.23.0.config``

.. code-block:: console

    # Local file list (directories)
    /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow/repo_code_v1.4_dbtest_0414.2025/dockerconfigs/fusion_multi_localdocker.v9.20241125.v0.23.0_v1.4.MANE.transcriptome.v3-1.config
    /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow/repo_code_v1.4_dbtest_0414.2025/dockerconfigs/fusion_multi_localdocker.v9.20241125.v0.23.0.config

----

----------------------
Area for improvements
----------------------

DB-update (Remaining tasks)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Task
     - Concept
     - Implementation

   * - V0.95 Probe annotation
     - Map 1039 probe regions to v0.95 pseudo transcriptome (Fusion V4 transcriptome)
     - Modify Steps 8-0-a to 8-0-c (loci file annotation), follow the “Fusion db_prep.steps_db-v3.1.docx” (ABIE-1012)

   * - V1.4 Kinase file check (update DB-v3.1 to DB-v3.2)
     - Check 26 MANE Select v1.4 transcripts to protein sequences mapping  
       
       Bug: BRAF in db-v3.1 uses "MANE v1.4 Plus clinical" (`ENST00000644969.2`) instead of correct "MANE v1.4 Select" (`ENST00000646891.2`)
     - Review the 26 kinase sequences (transcript ID, protein ID) manually curated in step 10 of `My PDF Guide <_static/Fusion_db_prep.steps_db-v3.1.pdf>`_

   * - V1.4 Whitelist update
     - Modify probe pair configuration:
       • AR:2,3,4 → AR:2,3  
       • BRAF:19-BRAF:11 → BRAF:18-BRAF:10
     - Review step 11 in `My PDF Guide <_static/Fusion_db_prep.steps_db-v3.1.pdf>`_

--------------------
Conclusion
--------------------

This completes the instructions for Fusion V4-based db construction.