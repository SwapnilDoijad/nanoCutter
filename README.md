# NanoCutter
NanoCutter is a tool to cut the Nanopore concatemeric and chimeric reads into individual repeat units.

Nanopore sequencing generates concatemeric or chimeric reads due to different reasons (circular DNA, erroneous ligation or failure of the pore/basecaller to split signals correctly).

## Prerequisites
- Python 3.9

## Installation
To install and set up NanoCutter, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/swapnildoijad/NanoCutter.git
   ```
2. Navigate to the project directory:
   ```bash
   cd NanoCutter
   ```
3. Make the scripts executable:
   ```bash
   chmod +x install.sh
   chmod +x ./NanoCutter.sh
   ```

4. Run the installation script:
   ```bash
   ./install.sh
   ```
5. Run test dataset:
   ```bash
   ./NanoCutter.sh -i test_data/fastq -o test_results_fastq
   ```

## Usage
```bash
Usage: ./NanoCutter.sh --help

Options:
        -h, --help               Show this help message and exit
        -c, --current-dir DIR    Set the current working directory (default: pwd)
        -i, --input-dir DIR      Set the input directory containing FASTA/FASTQ files
        -o, --output-dir DIR     Set the output directory for results

Examples:
        ./NanoCutter.sh -i test_data/fasta -o test_results_fasta
		./NanoCutter.sh -i test_data/fastq -o test_results_fastq
```
---
🧑‍💻 Author: Swapnil Doijad (swapnil.doijad@gmail.com)  
🙋 Support If you encounter bugs or have feature requests, please open an issue.