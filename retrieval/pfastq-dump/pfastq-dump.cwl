cwlVersion: v1.0
class: CommandLineTool
label: "pfastq-dump: A bash implementation of parallel-fastq-dump, parallel fastq-dump wrapper"
doc: "pfastq-dump is a bash implementation of parallel-fastq-dump, parallel fastq-dump wrapper. --stdout option is additionally supported, but almost same features. It also uses -N and -X options of fastq-dump to specify blocks of data to be decompressed separately. https://github.com/inutano/pfastq-dump"

$namespaces:
  edam: 'http://edamontology.org/'

hints:
  DockerRequirement:
    dockerPull: quay.io/inutano/sra-toolkit:v2.9.0

baseCommand: [pfastq-dump]

inputs:
  sraFiles:
    type: File[]
    inputBinding:
      position: 50
  nthreads:
    type: int
    inputBinding:
      prefix: -t
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
    format: edam:format_1930
    outputBinding:
      glob: "*fastq*"
  forward:
    type: File?
    format: edam:format_1930
    outputBinding:
      glob: "*_1.fastq*"
  reverse:
    type: File?
    format: edam:format_1930
    outputBinding:
      glob: "*_2.fastq*"
