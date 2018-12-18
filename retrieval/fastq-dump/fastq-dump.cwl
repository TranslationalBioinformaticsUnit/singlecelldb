cwlVersion: v1.0
class: CommandLineTool
label: "fastq-dump: dump .sra format file to generate fastq file"
doc: "sra-toolkit: https://github.com/ncbi/sra-tools/wiki/Download-On-Demand"

hints:
  DockerRequirement:
    dockerPull: quay.io/inutano/sra-toolkit:v2.9.0

baseCommand: [fastq-dump]

inputs:
  sraFiles:
    type: File[]
    inputBinding:
      position: 50
  split_files:
    type: boolean?
    default: true
    inputBinding:
      prefix: --split-files
  split_spot:
    type: boolean?
    default: true
    inputBinding:
      prefix: --split-spot
  skip_technical:
    type: boolean?
    default: true
    inputBinding:
      prefix: --skip-technical
  readids:
    type: boolean?
    default: true
    inputBinding:
      prefix: --readids
  gzip:
    type: boolean?
    default: true
    inputBinding:
      prefix: --gzip

outputs:
  fastqFiles:
    type: File[]
    outputBinding:
      glob: "*fastq*"
  forward:
    type: File?
    outputBinding:
      glob: "*_1.fastq*"
  reverse:
    type: File?
    outputBinding:
      glob: "*_2.fastq*"
