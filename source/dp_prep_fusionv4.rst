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

----------------------------
Executed commands (db-v3.1)
----------------------------

.. code-block:: console

    # !/bin/bash
    # code summary for db-v3.1 preparation

    ### Steps 1-2 download required files
    # link (mane 1.4): https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/
    # link (gencode): http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/
    # local_mane_dir="/mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/" # local MANE v1.4 db directory
    # local_gencode_dir="/mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/gencode_v47/" # local GENCODE v47 db directory

    ### Step 1
    # download MANE v1.4 source files 
    mkdir -p mane_v1.4 && cd mane_v1.4
    wget -e robots=off -r -np --no-check-certificate "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/"
    # move 1.4 files to another folder
    rsync -av ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4 OpenDB_MANE_human_v1.4
    rm -rf ftp.ncbi.nlm.nih.gov
    # mane v1.4 folder
    mkdir ./OpenDB_MANE_human_v1.4/derived/

    ### Step 2
    # derived files generation 
    ## namemap file
    zcat ./OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz > ./OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.summary.txt

    awk -F"\t" '{if($10 == "MANE Select")print $8"\t"$2"\t"$4"\t"$6}' ./OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.summary.txt > ./OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.summary.namemap

    ## gff file generation
    zcat ./mane_v1.4/OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gff.gz > ./mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.gff

    ## transcript list generation 
    mkdir -p ./Output_Final/
    awk '{print $1}' ./OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.summary.namemap > ./Output_Final/PA053_ACTFusionV5_PseudoIntron_MANE-v1.4_GENCODE-r47_capture-v1.0_GRCh38.20250407.transcript.MANE.list

    # download gencode v47 and move the files to a local folder “OpenDB_GENCODE_human_r47”
    wget -e robots=off -r -np http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/
    rsync -av ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/ OpenDB_GENCODE_human_r47
    rm -rf ftp.ebi.ac.uk
    mkdir ./OpenDB_GENCODE_human_r47/derived/

    ## generate derived files
    zcat ./OpenDB_GENCODE_human_r47/GRCh38.p14.genome.fa.gz > ./OpenDB_GENCODE_human_r47/derived/GRCh38.p14.genome.fa

    zcat ./OpenDB_GENCODE_human_r47/gencode.v47.annotation.gff3.gz > ./OpenDB_GENCODE_human_r47/derived/gencode.v47.annotation.gff3

    # actgenomics/fusion_dev:v0.6
    samtools faidx ./OpenDB_GENCODE_human_r47/derived/GRCh38.p14.genome.fa

    ### Step 3: retrieve transcript gff file
    python3 /mnt/RD_Develop/sandyteng/ACTFusionV5/code/filter_mane_gff.py \
    -i /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gff.gz \
    -o /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.transcript.gff

    ### Step 4: gff to bed conversion with ("bedops_2.4.39/bin/convert2bed")
    /tools/Fusion/convert2bed -i gff -d \
    < /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.gff > /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.bed

    bedtools getfasta -name -s \
    -fi ./gencode_v47/OpenDB_GENCODE_human_r47/derived/GRCh38.p14.genome.fa \
    -bed ./mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.bed \
    -fo ./mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.corrected.strand.fasta

    sed -i 's/([+-])//g' ./mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.corrected.strand.fasta

    ### Step 5:  generate plain annotation files
    python3 /mnt/RD_Develop/sandyteng/ACTFusionV5/code/RefFusion.v2.py \
    -g /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.gff \
    -m /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.summary.txt \
    -f /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.ensembl_genomic.transcript.corrected.strand.fasta \
    -p /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/Output_Final/PA053_ACTFusionV5_PseudoIntron_MANE-v1.4_GENCODE-r47_capturev1.0_GRCh38.20250407.transcript.MANE.only.list \
    -o /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47

    ### Step 8-0-a
    # obtain mapping exons (pseudo locations on 10*N transcriptome)
    bash /mnt/RD_Develop/sandyteng/ACTFusionV5/code/candidate_exons_mapping.sh \
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47.genome.loci \
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47.transcript.loci \
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.select.and.plus.clinical.namemap \
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/captureprobe_250401/ACTFusionv5_target-region_PartAB_individual_1039.bed \
    fusionv4.MANE.v1.4.GENCODE.r47 \
    /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/ \
    /tools/Fusion

    # extract mapped exons (candidate.exons.transcript.bed) sequences from gencode fasta file (gencode.genome.fa)
    bedtools getfasta -name -s \
    -fi /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47.fasta \
    -bed /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/fusionv4.MANE.v1.4.GENCODE.r47.candidate.exons.transcript.bed \
    -fo /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/probeseq/MANE.GRCh38.v1.4.0407.probe.r47.fasta

    ### Step 8-0-b
    # probe fasta generation
    python3 /mnt/RD_Develop/sandyteng/ACTFusionV5/code/Probe_faheader_converter.py \
    -f /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/probeseq/MANE.GRCh38.v1.4.0407.probe.r47.fasta \
    -n /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.select.and.plus.clinical.namemap \
    -o /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/probeseq/MANE.GRCh38.v1.4.0407.r47.probe.wtprimerlikeheader.fasta.gz

    # unzip fasta.gz
    gunzip /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/probeseq/MANE.GRCh38.v1.4.0407.r47.probe.wtprimerlikeheader.fasta.gz

    ### Step 8-0-c
    # perform probe sequences to transcriptome alignment
    /tools/Fusion/ncbi-blast/bin/blastn \
    -query /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/probeseq/MANE.GRCh38.v1.4.0407.r47.probe.wtprimerlikeheader.fasta \
    -subject /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47.fasta -outfmt 6 -task blastn-short -dust no \
    > /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/blastn/20250407_probe.r47.blastn

    # create blastn result for “reverse probe” and concatenate all the alignments
    cat 20250407_probe.r47.blastn > 20250407_rprobe.r47.blastn
    sed -i 's/|F|/|R|/' 20250407_rprobe.r47.blastn
    sed -i 's/mane/rmane/' 20250407_rprobe.r47.blastn
    cat 20250407_probe.r38.blastn 20250407_rprobe.r38.blastn > 20250407_probe.rprobe.r47.blastn

    # blastn parser (loci annotation)
    python3 /mnt/RD_Develop/sandyteng/ACTFusionV5/code/blastnparser.py \
    -if /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/InhouseDB_Probe/captureprobe_250407_MANE_Select/blastn/20250407_probe.rprobe.r47.blastn \
    -mp /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/mane_v1.4/OpenDB_MANE_human_v1.4/derived/MANE.GRCh38.v1.4.select.and.plus.clinical.namemap \
    -lf /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.r47.transcript.loci > /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_Loci/250407/PA053_ACTFusionV5_PseudoIntron_MANE-v1.4_GENCODE-r47_capturev1.0_GRCh38.20250407.transcript.MANE.only.blastn.r47.loci

    ### Step 9
    /tools/Fusion/bwa index /mnt/RD_Develop/sandyteng/ACTFusionV5/db_fusionv5/Output_MANE_Select/20250407_MANE.fasta

    ### Step 10 (manually curated from website)

    ### Step 11 ('filter_internal.QC9.0.mgsp.qcr.0.5-dbv3.v1.4.config' was modified in-place)
    cp /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/testconfigs/filter_internal.QC9.0.mgsp.qcr.0.5.blank.config /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/testconfigs/filter_internal.QC9.0.mgsp.qcr.0.5-dbv3.v1.4.config

    python3 /mnt/BI3/Team_workdir/sandyteng_workdir/ACTFusionV4_Torrent/code/update_qcconfig_with_tsv.py \
    -f /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/testconfigs/filter_internal.QC9.0.mgsp.qcr.0.5-dbv3.v1.4.config \
    -t /mnt/RD_Develop/sandyteng/ACTFusionV5/test/20250423_fusionv42v5_whitelist_gsppair/data/gsppairs_inclusion_v1.4.txt

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