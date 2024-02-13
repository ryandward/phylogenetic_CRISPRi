# Phylogenetic CRISPRi

This repository is dedicated to the development of CRISPRi systems for enteric gram-negative bacteria. It aims to streamline the integration of CRISPRi into the genetic makeup of these organisms, facilitating advanced research and applications.

## Detailed Information
For project documentation and more detailed information, please visit our [wiki page](https://github.com/ryandward/phylogenetic_CRISPRi/wiki).

## Custom Tools
The tools developed for this project, including sequence processing and analysis utilities, are hosted at [barcoder](https://github.com/ryandward/barcoder). Among these tools is `heuristicount`, a component of `barcoder` designed for accurate read counting. It analyzes the orientation, size, position, and valid context of reads, efficiently processing both `.fastq` and `.reads.zst` files.

## Data Storage
Sequences within this repository are compressed using a custom format, `*.reads.zst`. This method, implemented via `barcoder.distillreads`, not only significantly reduces file sizes to 5-10% of their original volume but also applies pairwise sorting of reads to further reduce their size. This format encapsulates a streamlined list of reads in `.zst` compression, omitting the additional metadata typically found in `.fastq` files.

## Contributions and Inquiries
For further inquiries or contributions, please refer to the [wiki](https://github.com/ryandward/phylogenetic_CRISPRi/wiki) or raise an issue in the repository.
