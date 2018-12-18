cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:

  fq: File
  genome_index: File
  down_trans_assem: boolean
  enable_cufflinks_output: boolean 

  hisat2_result_file: string

  bs_option: boolean
  samtools-view_result_file: string

  samtools-sort_result_file: string

  samtools-index_result_file: string

  annotation_file: File

  stringtie_result_file: string

  thread: int?

outputs:
  hisat2_version_stdout_result:
    type: File
    outputSource: hisat2_stdout/hisat2_version_stdout
  hisat2_version_result:
    type: File
    outputSource: hisat2_version/version_output
  hisat2_result:
    type: File
    outputSource: hisat2/hisat2_sam
  samtools_version_stderr_result:
    type: File
    outputSource: samtools_stderr/samtools_version_stderr
  samtools_version_result:
    type: File
    outputSource: samtools_version/version_output
  samtools-view_result:
    type: File
    outputSource: samtools-view/samtools-view_bam
  samtools-sort_result:
    type: File
    outputSource: samtools-sort/samtools-sort_sortbam
  samtools-index_result:
    type: File
    outputSource: samtools-index/samtools-index_indexbam
  cufflinks_version_stderr_result:
    type: File
    outputSource: cufflinks_stderr/cufflinks_version_stderr
  cufflinks_version_result:
    type: File
    outputSource: cufflinks_version/version_output
  cufflinks_result:
    type:
      type: array
      items: File
    outputSource: cufflinks/cufflinks_result
  stringtie_version_stdout_result:
    type: File
    outputSource: stringtie_stdout/stringtie_version_stdout
  stringtie_version_result:
    type: File
    outputSource: stringtie_version/version_output
  stringtie_result:
    type: File
    outputSource: stringtie/stringtie_result

steps:
  hisat2_stdout:
    run: hisat2-version.cwl
    in: []
    out: [hisat2_version_stdout]
  hisat2_version:
    run: ngs-version.cwl
    in:
      infile: hisat2_stdout/hisat2_version_stdout
    out: [version_output]
  hisat2:
    run: hisat2-se.cwl
    in:
      fq: fq
      gi: genome_index
      dta: down_trans_assem
      dta_cufflinks: enable_cufflinks_output
      output: hisat2_result_file
      process: thread
    out: [hisat2_sam]
  samtools_stderr:
    run: samtools-version.cwl
    in: []
    out: [samtools_version_stderr]
  samtools_version:
    run: ngs-version.cwl
    in:
      infile: samtools_stderr/samtools_version_stderr
    out: [version_output]
  samtools-view:
    run: samtools-view.cwl
    in:
      bs: bs_option
      bam: samtools-view_result_file
      sam: hisat2/hisat2_sam
    out: [samtools-view_bam]
  samtools-sort:
    run: samtools-sort.cwl
    in:
      bam: samtools-view/samtools-view_bam
      sortbam: samtools-sort_result_file
    out: [samtools-sort_sortbam]
  samtools-index:
    run: samtools-index.cwl
    in:
      sortbam: samtools-sort/samtools-sort_sortbam
      indexbam: samtools-index_result_file
    out: [samtools-index_indexbam]
  cufflinks_stderr:
    run: cufflinks-version.cwl
    in: []
    out: [cufflinks_version_stderr]
  cufflinks_version:
    run: ngs-version.cwl
    in:
      infile: cufflinks_stderr/cufflinks_version_stderr
    out: [version_output]
  cufflinks:
    run: cufflinks.cwl
    in:
      annotation: annotation_file
      bam: samtools-sort/samtools-sort_sortbam
      process: thread
    out: [cufflinks_result]
  stringtie_stdout:
    run: stringtie-version.cwl
    in: []
    out: [stringtie_version_stdout]
  stringtie_version:
    run: ngs-version.cwl
    in:
      infile: stringtie_stdout/stringtie_version_stdout
    out: [version_output]
  stringtie:
    run: stringtie.cwl
    in:
      bam: samtools-sort/samtools-sort_sortbam
      annotation: annotation_file
      output: stringtie_result_file
    out: [stringtie_result]
