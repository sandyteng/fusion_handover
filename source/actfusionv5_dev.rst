=================================
ACTFusion V5 Development Summary
=================================

-----------------
Purpose
-----------------

This document records the current status of ACTFusion V5 development. (data analysis summary and details)

----

-----------------
Relevant issues
-----------------

- `[Fusion V5] ACTFusion V5 Data Analysis summary (Sandy's handover) <https://actg.atlassian.net/browse/ABIE-1033>`_

----

------------------------------------------------
Executed (test/db construction/sequencing) runs
------------------------------------------------

ACTFusion V5 workflow execution summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: ACTFusion Workflow Test Scenarios
   :widths: 5 15 10 20 50
   :header-rows: 1

   * - Scenario
     - Workflow
     - Execution date
     - Type (Analysis)
     - Test summary
   * - 0 
     - Arriba workflow 
     - 2025/3/26 
     - Workflow construction test 
     - Sample for pipeline construction. See issue: ABIE-996
   * - 1 
     - Arriba workflow 
     - 2025/3/26 
     - Workflow construction test 
     - Most of the in-silico (IVT) fusions are not reported by Arriba caller.
   * - 2 
     - Arriba workflow 
     - 2025/3/26 
     - Workflow construction test 
     - The results (fusions.tsv vs fusions.tsv) are identical to /mnt/RD_Develop/sandyteng/ACTFusionV5/20250122_TwistBioscience/testresult/arriba_grch38/. Executed via Arriba’s containerized workflow.
   * - 3 
     - Fusion v4 workflow 
     - 2025/3/26 
     - Workflow construction test 
     - The workflow is copied from ACTFusion V4 (pipeline version: v0.29.0).
   * - 4 
     - Arriba workflow 
     - 2025/3/28 
     - Arriba test (300x IVT RNA) 
     - 300x IVT RNA test
   * - 5 
     - Fusion v4 workflow 
     - 2025/3/28 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v0.95 db)
   * - 6 
     - Fusion v4 workflow 
     - 2025/4/14 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v1.4 db - config test)
   * - 7 
     - Fusion v4 workflow 
     - 2025/4/16 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v1.4 db - v1 config test)
   * - 8 
     - Fusion v4 workflow 
     - 2025/4/21 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v1.4 db - v2 config test)
   * - 9 
     - Fusion v4 workflow 
     - 2025/4/21 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v1.4 db - v3 config test)
   * - 10 
     - Fusion v4 workflow 
     - 2025/4/23 
     - White list building (for db-v3.1) 
     - syn_exonskipping_seq-2-exons (63-9 = 53 reads) (MANE v1.4 db - v3 config test)
   * - 11 
     - Fusion v4 workflow 
     - 2025/4/23 
     - Fusion v1.4 db construction test 
     - 300x IVT RNA test (MANE v1.4 db - v3-1 config test)
   * - 12 
     - Fusion v4 workflow 
     - 2025/4/28 
     - Sequencing run analysis 
     - AANB01_504 (16 samples) (MANE v0.95 db)
   * - 13 
     - Arriba workflow 
     - 2025/4/30 
     - Sequencing run analysis 
     - AANB01_504 (16 samples) (GRCh38+RefSeq)
   * - 14 
     - Fusion v4 workflow 
     - 2025/5/5 
     - Sequencing run analysis 
     - AANB01_504 (16 samples) (MANE v1.4 db, db-v3.1)
   * - 15 
     - Fusion v4 workflow 
     - 2025/5/13 
     - Sequencing run analysis 
     - AANB01_507 (8 samples) (MANE v1.4 db, db-v3.1)
   * - 16 
     - Fusion v4 workflow 
     - 2025/5/13 
     - Sequencing run analysis 
     - AANB01_507 (8 samples) (MANE v0.95 db)
   * - 17 
     - Arriba workflow 
     - 2025/5/13 
     - Sequencing run analysis 
     - AANB01_507 (8 samples) (GRCh38+RefSeq)

Local nextflow output directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: console

    # Local directories
    /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow_outdir/
    ├── 20250325_PA043ANA_IVTALL-1_arriba_test
    ├── 20250326_PA043ANA_IVTRNA_arriba_test
    ├── 20250326_PA043ANA_IVTRNA_fusionv4_test
    ├── 20250326_Twist_8_NextSeq_samples_arriba_test
    ├── 20250328_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_arriba_test
    ├── 20250328_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_test
    ├── 20250414_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_MANE_v1.4_config_test
    ├── 20250416_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_MANE_v1.4_config_testv1
    ├── 20250418_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_MANE_v1.4_config_testv2
    ├── 20250421_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_MANE_v1.4_config_testv3
    ├── 20250423_ACTFusion_v4_81-IVTRNA-240906-300x_v0.1_fusionv4_MANE_v1.4_config_testv3-1
    ├── 20250423_syn_exonskipping_seq-2-exons_v0.1_fusionv4_MANE_v1.4_config_testv3
    ├── 20250428_AANB01_504-fusionv4-based_16_sample
    ├── 20250430_AANB01_504-arriba-based_16_sample
    ├── 20250505_AANB01_504-fusionv4-db-v3.1-based_16_sample
    ├── 20250512_AANB01_507-arriba-based_8_sample
    ├── 20250512_AANB01_507-fusionv4-based_8_sample
    └── 20250512_AANB01_507-fusionv4-db-v3.1-based_8_sample

    18 directories

--------------------
Conclusion
--------------------

This completes the record for ACTFusion V5 development summary.