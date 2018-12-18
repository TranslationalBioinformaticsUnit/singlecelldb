cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:

  run_id: string
  read_type: string

  genome_directory: Directory
  gtf: File
  output_sam_type: string
  output_name_prefix: string
  output_sam_strand_field: string

  annotation_file: File

  stringtie_result_file: string

  thread: int?

outputs:
  pfastq-dump_fq1_result:
    type: File
    outputSource: pfastq-dump/dump_fq1
  pfastq-dump_fq2_result:
    type: File
    outputSource: pfastq-dump/dump_fq2
  pfastq-dump_version_result:
    type: File
    outputSource: pfastq-dump/pfastq-dump_version
  star_version_stdout_result:
    type: File
    outputSource: star_stdout/star_version_stdout
  star_version_result:
    type: File
    outputSource: star_version/version_output
  star_result:
    type: File
    outputSource: star/star_bam
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
  pfastq-dump:
    run: prefetch_pfastq-dump.cwl
    in:
      run_id: run_id
      read_type: read_type
      thread: thread
    out: [dump_fq1, dump_fq2, pfastq-dump_version]
  star_stdout:
    run: star-version.cwl
    in: []
    out: [star_version_stdout]
  star_version:
    run: ngs-version.cwl
    in:
      infile: star_stdout/star_version_stdout
    out: [version_output]
  star:
    run: star-pe.cwl
    in:
      fq1: pfastq-dump/dump_fq1
      fq2: pfastq-dump/dump_fq2
      genome_directory: genome_directory
      gtf: gtf
      out_sam_type: output_sam_type
      out_name_prefix: output_name_prefix
      out_sam_strand_field: output_sam_strand_field
      process: thread
    out: [star_bam]
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
      bam: star/star_bam
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
      bam: star/star_bam
      annotation: annotation_file
      output: stringtie_result_file
    out: [stringtie_result]
