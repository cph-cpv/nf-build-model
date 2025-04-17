# `nf-build-reference`

_WORK IN PROGRESS_

A Nextflow workflow for building Iimi machine learning models.

1. `collapse_reference`

   Drop highly similar sequences from the build.

3. `find_unreliable_regions`

   Find regions of the reference sequences that are similar to a host genome.

4. `map_samples`

   Map a collection of sample FASTQ files against a reference FASTA.

## Subworkflows

### `collapse_reference`

Use `cd-hit-est` to cluster 99% similar sequences and produce a FASTA file of cluster representatives for sample mapping.

#### Outputs

* `reps.fa`

  The cluster sequences of the cluster representatives. The FASTA headers are UUIDs identifying the represenatives.

* `sequences_to_reps.csv`

  A CSV that maps sequence IDs (`sequence_id`) to their representatives (`rep_id`).

  ```text
  sequence_id,rep_id
  7e8ecf70-efb4-4adb-8b68-f54454785822,18582545-74f1-4e0e-be1c-2c964513f077
  59d6295c-18fb-4fa6-82f1-bca1666ae919,8873f6f4-3d86-470f-8532-0fd4c248932a
  f4d06f98-9a5d-4a66-8d96-e2fef66bac24,8873f6f4-3d86-470f-8532-0fd4c248932a
  ```

* `summary.txt`

  A summary of the changes to the reference sequences.

### `find_unreliable_regions`

Create pseudoreads from a reference FASTA file using a sliding window of 75.

Map pseudoreads against a host genome.

#### Inputs

* `host_fasta`

  A host FASTA file to map reference fragments against. This file will be used to build a Bowtie2 index.

* `reference_fasta`

  A FASTA file containing sequences to fragment and map against the host.

  The sequences in the file will be used to generated 75-mer fragments using a sliding window.

#### Outputs

* `unreliable_regions.csv`

  A CSV listing unreliable regions in all sequences.
  
  ```text
  sequence_id,start,end
  abcd1234,112,124
  abcd1234,598,612
  wxyz9876,78,124
  ```

### `map_samples`

#### Inputs

* `fasta`

  The FASTA file map against. A Bowtie2 index will be built from this file.

* `samples`

  A channel of sample FASTQ paths.

  The stem of each file is assumed to be the sample name. For example, `2023_SP_12345.fq.gz` will have the sample name `2023_SP_12345`.
