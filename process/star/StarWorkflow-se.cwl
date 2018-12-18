cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:

  fq: File
  genome_directory: Directory
  gtf: File
  output_sam_type: string
  output_name_prefix: string
  output_sam_strand_field: string

  annotation_file: File

  stringtie_result_file: string

  thread: int?

outputs:
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
    run: star-se.cwl
    in:
      fq: fq
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
