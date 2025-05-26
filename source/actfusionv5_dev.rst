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

--------------------------------
Executed (test/sequencing) runs
--------------------------------

Test Scenarios
==============

.. list-table:: Test scenarios for fusion workflows
   :header-rows: 1
   :widths: 8 20 12 40 60

   * - Scenario
     - Workflow
     - Execution date
     - Result directory
     - Test summary

   * - 0
     - Arriba workflow
     - 2025/03/26
     - /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow_outdir/20250325_PA043ANA_IVTALL-1_arriba_test/
     - Sample for pipeline construction. See issue: `ABIE-996 <https://actg.atlassian.net/browse/ABIE-996>`_

   * - 1
     - Arriba workflow
     - 2025/03/26
     - /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow_outdir/20250326_PA043ANA_IVTRNA_arriba_test/
     - Most of the in-silico (IVT) fusions are not reported by Arriba caller.

   * - 2
     - Arriba workflow
     - 2025/03/26
     - /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow_outdir/20250326_Twist_8_NextSeq_samples_arriba_test/
     - The results (``fusions.tsv`` vs ``fusions.tsv``) are identical to /mnt/RD_Develop/sandyteng/ACTFusionV5/20250122_TwistBioscience/testresult/arriba_grch38/. Executed via Arriba’s containerized workflow.

   * - 3
     - Fusion v4 workflow
     - 2025/03/26
     - /mnt/RD_Develop/sandyteng/ACTFusionV5/nextflow_outdir/20250326_PA043ANA_IVTRNA_fusionv4_test/
     - The workflow is copied from `ACTFusion V4 (pipeline version: v0.29.0) <https://github.com/ACTGenomics/torrent_fusion_pipeline_nextflow>`_.

--------------------
Conclusion
--------------------

This completes the record for ACTFusion V5 development summary.